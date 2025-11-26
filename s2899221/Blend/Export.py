import os
import json
import shutil
from mathutils import Vector
import bpy


# -----------------------------
# Configuration
# -----------------------------
# Export goes into the repo tree:
#   <home>/RayTracer/s2899221/ASCII/scene.json
#   <home>/RayTracer/s2899221/Textures/<images>
EXPORT_ROOT = os.path.join(os.path.expanduser("~"), "RayTracer", "s2899221")
TEXTURES_DIR = os.path.join(EXPORT_ROOT, "Textures")
SCENE_JSON_PATH = os.path.join(EXPORT_ROOT, "ASCII", "scene.json")


# Cache for image textures: maps a unique image key -> relative path ("Textures/<file>")
_IMAGE_CACHE = {}

# Global object counter to provide numeric IDs across all exported objects
_NEXT_OBJECT_ID = 0


def next_object_id():
    global _NEXT_OBJECT_ID
    oid = _NEXT_OBJECT_ID
    _NEXT_OBJECT_ID += 1
    return oid


# -----------------------------
# Utilities
# -----------------------------
def ensure_directories():
    os.makedirs(EXPORT_ROOT, exist_ok=True)
    os.makedirs(TEXTURES_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(SCENE_JSON_PATH), exist_ok=True)


def vec3(v):
    return [float(v[0]), float(v[1]), float(v[2])]


def get_final_resolution(scene):
    sx = int(scene.render.resolution_x * scene.render.resolution_percentage / 100)
    sy = int(scene.render.resolution_y * scene.render.resolution_percentage / 100)
    return [sx, sy]


def classify_mesh_object(obj):
    """
    Classify the mesh as one of: Cube, Sphere, Cylinder, Plane, or None if unsupported.
    Heuristic: Based on object/data names (common Blender primitive defaults).
    """
    name_l = obj.name.lower()
    data_l = obj.data.name.lower() if obj.data else ""
    source = f"{name_l} {data_l}"
    if "cube" in source:
        return "Cube"
    if "sphere" in source or "uvsphere" in source or "ico" in source:
        return "Sphere"
    if "cylinder" in source:
        return "Cylinder"
    if "plane" in source:
        return "Plane"
    return None


def get_world_corners_for_plane(obj):
    """
    Returns the world-space corners (as 4 vec3 lists) for a plane-like mesh object.
    Assumes the mesh has at least 4 vertices (default Plane has exactly 4).
    """
    corners = []
    mw = obj.matrix_world
    mesh = obj.data
    if mesh is None or len(mesh.vertices) < 4:
        return corners
    # Use the first 4 vertices in their current order
    for i in range(4):
        co = mw @ mesh.vertices[i].co
        corners.append(vec3(co))
    return corners


def material_principled_node(mat):
    if not mat or not mat.use_nodes or not mat.node_tree:
        return None
    for node in mat.node_tree.nodes:
        if isinstance(node, bpy.types.ShaderNodeBsdfPrincipled):
            return node
    return None


def find_upstream_color_value(input_socket, depth=0, max_depth=16, visited=None):
    """
    Try to evaluate a constant color upstream of a socket.
    Supports simple cases: ShaderNodeRGB, MixRGB with constant inputs.
    Returns (type, payload):
      ("RGB", [r,g,b]) or ("ImageTexture", "Textures/<path>") or (None, None) if unknown.
    """
    if visited is None:
        visited = set()
    if depth > max_depth or not input_socket:
        return (None, None)
    if input_socket.is_linked:
        for link in input_socket.links:
            node = link.from_node
            if node is None or node in visited:
                continue
            visited.add(node)
            # Image texture
            if isinstance(node, bpy.types.ShaderNodeTexImage) and node.image:
                rel_path = copy_or_save_image_to_textures(node.image)
                if rel_path:
                    return ("ImageTexture", rel_path)
            # Plain RGB node
            if isinstance(node, bpy.types.ShaderNodeRGB):
                col = getattr(node.outputs[0], "default_value", None)
                if col is not None and len(col) >= 3:
                    return ("RGB", [float(col[0]), float(col[1]), float(col[2])])
            # MixRGB with constant inputs
            if isinstance(node, bpy.types.ShaderNodeMixRGB):
                fac = getattr(node.inputs[0], "default_value", 0.5)
                def col_from_input(inp):
                    if inp.is_linked:
                        t, p = find_upstream_color_value(inp, depth + 1, max_depth, visited)
                        if t == "RGB":
                            return p
                        return None
                    dv = getattr(inp, "default_value", None)
                    if dv is not None and len(dv) >= 3:
                        return [float(dv[0]), float(dv[1]), float(dv[2])]
                    return None
                c1 = col_from_input(node.inputs[1])
                c2 = col_from_input(node.inputs[2])
                if c1 and c2:
                    f = float(fac)
                    r = (1.0 - f) * c1[0] + f * c2[0]
                    g = (1.0 - f) * c1[1] + f * c2[1]
                    b = (1.0 - f) * c1[2] + f * c2[2]
                    return ("RGB", [r, g, b])
            # Recurse through all inputs
            for inp in getattr(node, "inputs", []):
                t, p = find_upstream_color_value(inp, depth + 1, max_depth, visited)
                if t is not None:
                    return (t, p)
    return (None, None)


def find_upstream_image_tex_node(input_socket, depth=0, max_depth=16, visited=None):
    """
    Walk upstream from a socket and try to find a ShaderNodeTexImage.
    Returns the first found image texture node or None.
    """
    if visited is None:
        visited = set()
    if depth > max_depth:
        return None
    if not input_socket or not input_socket.is_linked:
        return None
    for link in input_socket.links:
        node = link.from_node
        if node is None or node in visited:
            continue
        visited.add(node)
        if isinstance(node, bpy.types.ShaderNodeTexImage):
            return node
        # Recurse through all inputs of this node
        for inp in getattr(node, "inputs", []):
            found = find_upstream_image_tex_node(inp, depth + 1, max_depth, visited)
            if found:
                return found
    return None


def copy_or_save_image_to_textures(image):
    """
    Ensure the given Blender image is written/copied into TEXTURES_DIR.
    Uses a cache to avoid repeated work. Returns the relative path "Textures/<file>".
    """
    if image is None:
        return None

    # Unique key per image – prefer absolute filepath if present, otherwise name
    key = None
    abs_src = None
    if image.source == 'FILE' and image.filepath:
        abs_src = bpy.path.abspath(image.filepath)
        key = os.path.normpath(abs_src)
    else:
        key = f"__packed__::{image.name_full}"

    if key in _IMAGE_CACHE:
        return _IMAGE_CACHE[key]

    # Determine destination filename
    ext = os.path.splitext(image.filepath if image.filepath else image.name_full)[1]
    if not ext:
        ext = ".png"
    safe_name = image.name_full.replace(os.sep, "_").replace(":", "_")
    dest_filename = f"{safe_name}{ext}"
    dest_abs = os.path.join(TEXTURES_DIR, dest_filename)
    dest_rel = os.path.join("Textures", dest_filename)

    # Write/copy
    try:
        if abs_src and os.path.exists(abs_src):
            if not os.path.exists(dest_abs):
                shutil.copy2(abs_src, dest_abs)
        else:
            # Packed or generated image: save to disk
            # Preserve original path values to avoid mutating user's file state
            orig_filepath_raw = image.filepath_raw
            try:
                # Use PNG as a safe default for packed/generated sources
                if not dest_abs.lower().endswith(".png"):
                    dest_abs = os.path.splitext(dest_abs)[0] + ".png"
                    dest_rel = os.path.splitext(dest_rel)[0] + ".png"
                image.filepath_raw = dest_abs
                image.file_format = 'PNG'
                image.save()
            finally:
                image.filepath_raw = orig_filepath_raw
    except Exception as e:
        print(f"[export] Failed to export image '{image.name_full}': {e}")
        return None

    _IMAGE_CACHE[key] = dest_rel
    return dest_rel


def extract_material_info(obj):
    """
    Extract material properties for the object, assuming Principled BSDF.
    Returns:
      {
        "baseColor": { "type": "RGB", "value": [r,g,b] } | { "type": "ImageTexture", "path": "Textures/..." },
        "alpha": float,
        "metallic": float,
        "roughness": float,
        "ior": float
      }
    or None if no material.
    """
    mat = None
    if obj.material_slots:
        for slot in obj.material_slots:
            if slot.material:
                mat = slot.material
                break
    if mat is None:
        return None

    principled = material_principled_node(mat)
    base_color_info = None
    alpha = 1.0
    metallic = 0.0
    roughness = 0.5
    ior = 1.45

    if principled:
        base_input = principled.inputs.get("Base Color")
        if base_input:
            # First try to evaluate a color or texture upstream
            t, payload = find_upstream_color_value(base_input)
            if t == "ImageTexture":
                base_color_info = {"type": "ImageTexture", "path": payload}
            elif t == "RGB":
                base_color_info = {"type": "RGB", "value": payload}
            else:
                # Fallback: direct image lookup
                img_node = find_upstream_image_tex_node(base_input)
                if img_node and img_node.image:
                    rel_path = copy_or_save_image_to_textures(img_node.image)
                    if rel_path:
                        base_color_info = {"type": "ImageTexture", "path": rel_path}
        if base_color_info is None and base_input and hasattr(base_input, "default_value") and base_input.default_value:
            dv = base_input.default_value
            base_color_info = {"type": "RGB", "value": [float(dv[0]), float(dv[1]), float(dv[2])]}
            if len(dv) >= 4:
                alpha = float(dv[3])

        alpha_input = principled.inputs.get("Alpha")
        if alpha_input and hasattr(alpha_input, "default_value") and alpha_input.default_value is not None:
            alpha = float(alpha_input.default_value)

        met_input = principled.inputs.get("Metallic")
        if met_input and hasattr(met_input, "default_value") and met_input.default_value is not None:
            metallic = float(met_input.default_value)

        rough_input = principled.inputs.get("Roughness")
        if rough_input and hasattr(rough_input, "default_value") and rough_input.default_value is not None:
            roughness = float(rough_input.default_value)

        ior_input = principled.inputs.get("IOR")
        if ior_input and hasattr(ior_input, "default_value") and ior_input.default_value is not None:
            ior = float(ior_input.default_value)
    else:
        # Fallback to material diffuse color if nodes unavailable
        if hasattr(mat, "diffuse_color") and mat.diffuse_color:
            col = mat.diffuse_color
            base_color_info = {"type": "RGB", "value": [float(col[0]), float(col[1]), float(col[2])]}
            if len(col) >= 4:
                alpha = float(col[3])

    if base_color_info is None:
        # No color info found; default to white
        base_color_info = {"type": "RGB", "value": [1.0, 1.0, 1.0]}

    return {
        "baseColor": base_color_info,
        "alpha": float(alpha),
        "metallic": float(metallic),
        "roughness": float(roughness),
        "ior": float(ior)
    }


# -----------------------------
# Exporters
# -----------------------------
def export_cameras(scene):
    out = {"Perspective": []}
    for obj in bpy.data.objects:
        if obj.type != 'CAMERA':
            continue
        cam = obj.data
        if cam.type != 'PERSP':
            continue

        oid = next_object_id()
        loc = vec3(obj.location)
        # Blender camera looks along -Z in local space
        forward = (obj.matrix_world.to_quaternion() @ Vector((0.0, 0.0, -1.0))).normalized()
        direction = vec3(forward)

        focal_length_mm = float(cam.lens)
        sensor_w_mm = float(cam.sensor_width)
        sensor_h_mm = float(cam.sensor_height)
        resolution = get_final_resolution(scene)

        out["Perspective"].append({
            "id": oid,
            "location": loc,
            "direction": direction,
            "focalLengthMM": focal_length_mm,
            "sensorWidthMM": sensor_w_mm,
            "sensorHeightMM": sensor_h_mm,
            "filmResolution": resolution
        })
    return out


def export_lights():
    out = {"Point": []}
    for obj in bpy.data.objects:
        if obj.type != 'LIGHT':
            continue
        light = obj.data
        if light.type != 'POINT':
            continue

        oid = next_object_id()
        out["Point"].append({
            "id": oid,
            "location": vec3(obj.location),
            # Store Blender energy as radiant intensity proxy
            "radiantIntensity": float(getattr(light, "energy", 0.0))
        })
    return out


def export_meshes():
    out = {
        "Cube": [],
        "Sphere": [],
        "Cylinder": [],
        "Plane": []
    }
    for obj in bpy.data.objects:
        if obj.type != 'MESH':
            continue
        kind = classify_mesh_object(obj)
        if not kind:
            continue

        oid = next_object_id()
        mat_info = extract_material_info(obj)

        if kind == "Plane":
            corners = get_world_corners_for_plane(obj)
            out["Plane"].append({
                "id": oid,
                "corners": corners,
                "material": mat_info
            })
        else:
            # TRS in world space (derived from matrix_world to capture parenting/constraints and correct rotation)
            mw = obj.matrix_world
            loc_v = mw.to_translation()
            rot_q = mw.to_quaternion()  # (w, x, y, z)
            scl_v = mw.to_scale()
            # Also provide Euler in a fixed order derived from world matrix for compatibility
            rot_e = mw.to_euler('XYZ')
            out[kind].append({
                "id": oid,
                "translation": [float(loc_v.x), float(loc_v.y), float(loc_v.z)],
                "rotationQuat": [float(rot_q.w), float(rot_q.x), float(rot_q.y), float(rot_q.z)],
                "rotationEuler": [float(rot_e.x), float(rot_e.y), float(rot_e.z)],  # radians, XYZ order from world matrix
                "scale": [float(scl_v.x), float(scl_v.y), float(scl_v.z)],
                "material": mat_info
            })
    return out


def assemble_scene_json():
    scene = bpy.context.scene
    data = {
        "Camera": export_cameras(scene),
        "Light": export_lights(),
        "Mesh": export_meshes()
    }
    return data


def write_scene_json(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def main():
    ensure_directories()
    data = assemble_scene_json()
    write_scene_json(SCENE_JSON_PATH, data)
    print(f"[export] Scene exported to: {SCENE_JSON_PATH}")


if __name__ == "__main__":
    main()


