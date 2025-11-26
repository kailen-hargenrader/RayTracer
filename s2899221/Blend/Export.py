import warnings
import math
import json
import os
try:
    import bpy
except ImportError:
    raise SystemExit("blender not found, must run script in blender")

class id_generator:
    def __init__(self):
        self.id = 0
    def get_id(self):
        self.id += 1
        return self.id

id_gen = id_generator()

def vec3_to_dict(v):
    return {"x": float(v[0]), "y": float(v[1]), "z": float(v[2])}

def list3_to_dict(v):
    return {"x": float(v[0]), "y": float(v[1]), "z": float(v[2])}

def vec2_to_dict(v):
    return {"x": int(v[0]), "y": int(v[1])}

def _ensure_dir(dir_path: str) -> None:
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path, exist_ok=True)

def _to_u8(x: float) -> int:
    if x is None:
        return 0
    if x < 0.0:
        x = 0.0
    if x > 1.0:
        x = 1.0
    return int(round(x * 255.0))

# ------- CACHES TO AVOID RE-DOING WORK -------
_material_params_cache = {}
_image_to_ppm_cache = {}

def _image_cache_key(img) -> str:
    # Prefer filepath key so duplicates of the same image datablock dedupe
    try:
        fp = getattr(img, "filepath_raw", None) or getattr(img, "filepath", None)
        if fp:
            return str(fp)
    except Exception:
        pass
    try:
        return str(getattr(img, "name", ""))
    except Exception:
        return ""

def save_image_as_ppm(img, out_path: str) -> str:
    """
    Save a Blender image datablock to binary PPM (P6).
    Returns the absolute output path on success, raises on failure.
    """
    # Fast-path: if we've already exported this image, reuse the cached path
    key = _image_cache_key(img)
    if key in _image_to_ppm_cache:
        return _image_to_ppm_cache[key]
    # Some images may be packed; ensure pixels are available
    w = int(getattr(img, "size", [0, 0])[0])
    h = int(getattr(img, "size", [0, 0])[1])
    if w <= 0 or h <= 0:
        raise RuntimeError("Image has invalid size")
    # Blender stores linear floats RGBA in pixels
    px = list(getattr(img, "pixels", []))
    if not px or len(px) < w * h * 4:
        # Attempt to load from disk if available
        try:
            img.update()
            px = list(getattr(img, "pixels", []))
        except Exception:
            pass
    if not px or len(px) < w * h * 4:
        raise RuntimeError("Could not access image pixel buffer")
    _ensure_dir(os.path.dirname(out_path))
    # If a file with the same name already exists, reuse it (assume unchanged)
    if os.path.isfile(out_path):
        abs_path = os.path.abspath(out_path)
        _image_to_ppm_cache[key] = abs_path
        return abs_path
    with open(out_path, "wb") as f:
        header = f"P6\n{w} {h}\n255\n".encode("ascii")
        f.write(header)
        # Write row-major top-to-bottom as stored
        # Convert float RGBA to 8-bit RGB (drop alpha)
        # Note: no gamma correction applied
        for i in range(0, w * h * 4, 4):
            r = _to_u8(px[i + 0])
            g = _to_u8(px[i + 1])
            b = _to_u8(px[i + 2])
            f.write(bytes((r, g, b)))
    abs_path = os.path.abspath(out_path)
    _image_to_ppm_cache[key] = abs_path
    return abs_path

def _export_albedo_texture_if_needed(img) -> str:
    """
    Export the given Blender image as PPM once and reuse the result.
    """
    img_name = getattr(img, "name", "texture")
    safe = "".join(c if c.isalnum() or c in ("_", "-", ".") else "_" for c in img_name)
    base_dir = os.path.join(os.path.expanduser("~"), "RayTracer", "s2899221", "Textures")
    out_ppm = os.path.join(base_dir, f"{safe}_basecolor.ppm")
    try:
        return save_image_as_ppm(img, out_ppm)
    except Exception:
        return ""

def extract_material_params_from_object(obj):
    alpha = 1.0
    metallic = 0.0
    roughness = 0.5
    ior = 1.5
    # Albedo defaults: neutral white (so legacy grayscale lighting stays visible)
    albedo = [1.0, 1.0, 1.0]
    albedo_texture_ppm = None

    try:
        # Prefer active material; otherwise first slot with a material
        material = getattr(obj, "active_material", None)
        if material is None:
            slots = getattr(obj, "material_slots", [])
            for slot in slots:
                m = getattr(slot, "material", None)
                if m is not None:
                    material = m
                    break

        # If this object has a material we have already processed exactly, reuse it
        if material is not None:
            try:
                mat_key = int(material.as_pointer())
            except Exception:
                mat_key = getattr(material, "name", None)
            if mat_key is not None and mat_key in _material_params_cache:
                return dict(_material_params_cache[mat_key])

        if material is not None:
            if getattr(material, "use_nodes", False) and getattr(material, "node_tree", None):
                nt = material.node_tree
                principled = None
                for node in getattr(nt, "nodes", []):
                    # Node.type is 'BSDF_PRINCIPLED' for Principled BSDF
                    if getattr(node, "type", "") == "BSDF_PRINCIPLED" or getattr(node, "bl_idname", "") == "ShaderNodeBsdfPrincipled":
                        principled = node
                        break

                if principled is not None:
                    def _get_input_value(node, name, fallback):
                        try:
                            sock = node.inputs.get(name)
                        except Exception:
                            sock = None
                        if sock is None:
                            return float(fallback)
                        try:
                            val = getattr(sock, "default_value", fallback)
                            if isinstance(val, (list, tuple)):
                                return float(val[0]) if len(val) > 0 else float(fallback)
                            return float(val)
                        except Exception:
                            return float(fallback)

                    # Metallic/roughness/alpha/IOR scalars
                    metallic = _get_input_value(principled, "Metallic", metallic)
                    roughness = _get_input_value(principled, "Roughness", roughness)
                    alpha = _get_input_value(principled, "Alpha", alpha)
                    ior = _get_input_value(principled, "IOR", ior)

                    # Base Color: prefer socket default if not linked
                    try:
                        base_color_socket = principled.inputs.get("Base Color")
                    except Exception:
                        base_color_socket = None
                    if base_color_socket is not None:
                        # If linked, try to find an Image Texture node and bake to PPM
                        links = list(getattr(base_color_socket, "links", []))
                        if links:
                            try:
                                from_node = links[0].from_node
                            except Exception:
                                from_node = None
                            if from_node and getattr(from_node, "bl_idname", "") == "ShaderNodeTexImage":
                                img = getattr(from_node, "image", None)
                                if img:
                                    # Export once per unique image and reuse
                                    exported = _export_albedo_texture_if_needed(img)
                                    if exported:
                                        albedo_texture_ppm = exported
                        # If not linked or baking failed, read default color (RGBA)
                        if albedo_texture_ppm is None:
                            try:
                                dv = getattr(base_color_socket, "default_value", [1.0, 1.0, 1.0, 1.0])
                                if isinstance(dv, (list, tuple)) and len(dv) >= 3:
                                    albedo = [float(dv[0]), float(dv[1]), float(dv[2])]
                            except Exception:
                                pass
            else:
                # Fallbacks for non-node materials if available
                try:
                    metallic = float(getattr(material, "metallic", metallic))
                except Exception:
                    pass
                try:
                    roughness = float(getattr(material, "roughness", roughness))
                except Exception:
                    pass
                try:
                    alpha = float(getattr(material, "alpha", alpha))
                except Exception:
                    pass
                try:
                    ior = float(getattr(material, "ior", ior))
                except Exception:
                    pass
                # Try to get diffuse_color as albedo if present
                try:
                    dc = getattr(material, "diffuse_color", None)
                    if isinstance(dc, (list, tuple)) and len(dc) >= 3:
                        albedo = [float(dc[0]), float(dc[1]), float(dc[2])]
                except Exception:
                    pass
    except Exception:
        # Silently keep defaults if anything goes wrong
        pass

    result = {
        "alpha": float(alpha),
        "metallic": float(metallic),
        "roughness": float(roughness),
        "ior": float(ior),
        # Store as vec3 in x,y,z so existing parser can read it
        "albedo": list3_to_dict(albedo),
        # Absolute path to baked PPM if available
        "albedo_texture_ppm": (albedo_texture_ppm or ""),
    }
    # Cache material params for reuse across objects sharing the same material
    try:
        if material is not None:
            mat_key = int(material.as_pointer())
            _material_params_cache[mat_key] = dict(result)
    except Exception:
        pass
    return result

def mesh_handler(mesh_dict, mesh_obj):
    def cube_handler(cube_dict, cube_obj):
        mw = cube_obj.matrix_world
        translation = vec3_to_dict(mw.to_translation())
        eul = mw.to_euler('XYZ')
        rotation = {
            "roll": float(eul.x),
            "pitch": float(eul.y),
            "yaw": float(eul.z),
        }
        scale_vec = mw.to_scale()
        scale = vec3_to_dict(scale_vec)

        material = extract_material_params_from_object(cube_obj)

        cube_params = {
            "translation": translation,
            "rotation": rotation,
            "scale": scale,
            "material": material,
        }
        cube_dict[id_gen.get_id()] = cube_params
    def cylinder_handler(cylinder_dict, cylinder_obj):
        mw = cylinder_obj.matrix_world
        translation = vec3_to_dict(mw.to_translation())
        eul = mw.to_euler('XYZ')
        rotation = {
            "roll": float(eul.x),
            "pitch": float(eul.y),
            "yaw": float(eul.z),
        }
        scale_vec = mw.to_scale()
        # Export full 3D scale vector only; consumers can derive radius/length if needed
        scale = vec3_to_dict(scale_vec)

        material = extract_material_params_from_object(cylinder_obj)

        cylinder_params = {
            "translation": translation,
            "rotation": rotation,
            "scale": scale,
            "material": material,
        }
        cylinder_dict[id_gen.get_id()] = cylinder_params
    def sphere_handler(sphere_dict, sphere_obj):
        mw = sphere_obj.matrix_world
        location = list(mw.to_translation())
        scale_vec = mw.to_scale()
        if abs(scale_vec.x - scale_vec.y) > 1e-6 or abs(scale_vec.y - scale_vec.z) > 1e-6:
            warnings.warn("Non-uniform scale on sphere; exporting non-uniform scale vector")

        material = extract_material_params_from_object(sphere_obj)

        sphere_params = {
            "location": vec3_to_dict(location),
            "scale": vec3_to_dict(scale_vec),
            "material": material,
        }
        sphere_dict[id_gen.get_id()] = sphere_params
    def plane_handler(plane_dict, plane_obj):
        mw = plane_obj.matrix_world
        me = plane_obj.data

        corners = []

        def _transform_point(m, p):
            x, y, z = p
            return [
                m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3],
                m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3],
                m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3],
            ]

        # Prefer the actual quad polygon corners if available
        for poly in getattr(me, "polygons", []):
            vert_indices = list(getattr(poly, "vertices", []))
            if len(vert_indices) == 4:
                corners = [_transform_point(mw, tuple(me.vertices[i].co)) for i in vert_indices]
                break

        # Fallback: derive 4 corners from the object's bounding box
        if not corners:
            bb_local = [tuple(c) for c in getattr(plane_obj, "bound_box", [])]
            if bb_local:
                bb_world = [_transform_point(mw, p) for p in bb_local]
                xs = [p[0] for p in bb_world]
                ys = [p[1] for p in bb_world]
                zs = [p[2] for p in bb_world]
                dx = max(xs) - min(xs)
                dy = max(ys) - min(ys)
                dz = max(zs) - min(zs)
                # Identify the thinnest axis and pick the corresponding face (4 points)
                eps = 1e-6
                if dx <= dy and dx <= dz:
                    minx, maxx = min(xs), max(xs)
                    candidates = [p for p in bb_world if abs(p[0] - minx) < eps]
                    if len(candidates) != 4:
                        candidates = [p for p in bb_world if abs(p[0] - maxx) < eps]
                elif dy <= dx and dy <= dz:
                    miny, maxy = min(ys), max(ys)
                    candidates = [p for p in bb_world if abs(p[1] - miny) < eps]
                    if len(candidates) != 4:
                        candidates = [p for p in bb_world if abs(p[1] - maxy) < eps]
                else:
                    minz, maxz = min(zs), max(zs)
                    candidates = [p for p in bb_world if abs(p[2] - minz) < eps]
                    if len(candidates) != 4:
                        candidates = [p for p in bb_world if abs(p[2] - maxz) < eps]
                if len(candidates) == 4:
                    corners = [list(p) for p in candidates]
                else:
                    warnings.warn("Could not determine plane corners from bounding box; exporting first 4 bbox points")
                    corners = [list(p) for p in bb_world[:4]]
            else:
                warnings.warn("Plane has no bounding box; corners unavailable")

        material = extract_material_params_from_object(plane_obj)

        plane_params = {
            "corners": [list3_to_dict(p) for p in corners],
            "material": material,
        }
        plane_dict[id_gen.get_id()] = plane_params

    SUPPORTED_MESH_DATATYPE_TO_HANDLER = {
        "Cube": cube_handler,
        "Cylinder": cylinder_handler,
        "Sphere": sphere_handler,
        "Plane": plane_handler,
    }

    def is_supported_mesh_type(mesh_type):
        return mesh_type in SUPPORTED_MESH_DATATYPE_TO_HANDLER.keys()

    mesh_type = mesh_obj.data.name
    mesh_type = mesh_type.split('.')[0] # remove .001, .002, etc.
    if not is_supported_mesh_type(mesh_type):
        warnings.warn(f"Unsupported mesh type, skipping object: {mesh_type}")
    else:
        if mesh_type not in mesh_dict.keys():
            mesh_dict[mesh_type] = {}
        type_dict = mesh_dict[mesh_type]
        mesh_handler_func = SUPPORTED_MESH_DATATYPE_TO_HANDLER[mesh_type]
        mesh_handler_func(type_dict, mesh_obj)


def camera_handler(cam_dict, cam_obj):
    def perspective_handler(persp_dict, persp_obj):
        mw = persp_obj.matrix_world
        # World-space camera location
        location = vec3_to_dict(mw.to_translation())
        # Camera looks down local -Z. Transform to world and normalize
        gaze = [-mw[0][2], -mw[1][2], -mw[2][2]]
        length = (gaze[0] ** 2 + gaze[1] ** 2 + gaze[2] ** 2) ** 0.5
        if length > 1e-8:
            direction = [gaze[0] / length, gaze[1] / length, gaze[2] / length]
        else:
            direction = [0.0, 0.0, -1.0]

        cam_data = getattr(persp_obj, "data", None)
        focal_length = float(getattr(cam_data, "lens", 0.0))
        sensor_width = float(getattr(cam_data, "sensor_width", 0.0))
        sensor_height = float(getattr(cam_data, "sensor_height", 0.0))

        # Film resolution from render settings (apply percentage scale)
        try:
            res_x = int(scene.render.resolution_x * scene.render.resolution_percentage / 100)
            res_y = int(scene.render.resolution_y * scene.render.resolution_percentage / 100)
        except Exception:
            res_x, res_y = 0, 0
        film_resolution = vec2_to_dict([res_x, res_y])

        persp_params = {
            "location": location,
            "direction": list3_to_dict(direction),
            "focal_length": focal_length,
            "sensor_width": sensor_width,
            "sensor_height": sensor_height,
            "film_resolution": film_resolution,
        }
        persp_dict[id_gen.get_id()] = persp_params

    SUPPORTED_CAMERA_DATATYPE_TO_HANDLER = {
        "PERSP": perspective_handler,
    }

    def is_supported_cam_type(cam_type):
        return cam_type in SUPPORTED_CAMERA_DATATYPE_TO_HANDLER.keys()

    cam_type = cam_obj.data.type
    if not is_supported_cam_type(cam_type):
        warnings.warn(f"Unsupported mesh type, skipping object: {cam_type}")
    else:
        if cam_type not in cam_dict.keys():
            cam_dict[cam_type] = {}
        type_dict = cam_dict[cam_type]
        cam_handler_func = SUPPORTED_CAMERA_DATATYPE_TO_HANDLER[cam_type]
        cam_handler_func(type_dict, cam_obj)

def light_handler(light_dict, light_obj):
    def point_handler(point_dict, point_obj):
        mw = point_obj.matrix_world
        location = vec3_to_dict(mw.to_translation())
        energy = float(getattr(point_obj.data, "energy", 0.0))
        radiant_intensity = energy / (4.0 * math.pi)

        point_params = {
            "location": location,
            "radiant_intensity": radiant_intensity,
        }
        point_dict[id_gen.get_id()] = point_params

    SUPPORTED_LIGHT_DATATYPE_TO_HANDLER = {
        "POINT": point_handler,
    }

    def is_supported_light_type(light_type):
        return light_type in SUPPORTED_LIGHT_DATATYPE_TO_HANDLER.keys()

    light_type = light_obj.data.type
    if not is_supported_light_type(light_type):
        warnings.warn(f"Unsupported mesh type, skipping object: {light_type}")
    else:
        if light_type not in light_dict.keys():
            light_dict[light_type] = {}
        type_dict = light_dict[light_type]
        light_handler_func = SUPPORTED_LIGHT_DATATYPE_TO_HANDLER[light_type]
        light_handler_func(type_dict, light_obj)


SUPPORTED_OBJECT_TYPE_TO_HANDLER = {
    "MESH": mesh_handler,
    "CAMERA": camera_handler,
    "LIGHT": light_handler,
}

scene = bpy.context.scene

def get_data_from_scene(scene):
    def is_supported_obj(obj_type):
        return obj_type in SUPPORTED_OBJECT_TYPE_TO_HANDLER.keys()

    data = {}
    for obj in scene.objects:
        obj_type = getattr(obj, "type", type(obj).__name__)
        if not is_supported_obj(obj_type):
            warnings.warn(f"Unsupported object type, skipping object: {obj_type}")
        else:
            if obj_type not in data.keys():
                data[obj_type] = {}
            obj_dict = data[obj_type]
            obj_handler_func = SUPPORTED_OBJECT_TYPE_TO_HANDLER[obj_type]
            obj_handler_func(obj_dict, obj)
    return data

def export_json_from_data(data, path):
    json_data = json.dumps(data, indent=4)
    with open(path, "w") as f:
        f.write(json_data)


if __name__ == "__main__":
    scene = bpy.context.scene
    data = get_data_from_scene(scene)
    filepath = os.path.join(os.path.expanduser("~"), "RayTracer", "s2899221", "ASCII", "cylinders.json")
    export_json_from_data(data, filepath)
