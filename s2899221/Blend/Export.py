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

        cube_params = {
            "translation": translation,
            "rotation": rotation,
            "scale": scale,
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

        cylinder_params = {
            "translation": translation,
            "rotation": rotation,
            "scale": scale,
        }
        cylinder_dict[id_gen.get_id()] = cylinder_params
    def sphere_handler(sphere_dict, sphere_obj):
        mw = sphere_obj.matrix_world
        location = list(mw.to_translation())
        scale_vec = mw.to_scale()
        if abs(scale_vec.x - scale_vec.y) > 1e-6 or abs(scale_vec.y - scale_vec.z) > 1e-6:
            warnings.warn("Non-uniform scale on sphere; exporting non-uniform scale vector")

        sphere_params = {
            "location": vec3_to_dict(location),
            "scale": vec3_to_dict(scale_vec),
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

        plane_params = {
            "corners": [list3_to_dict(p) for p in corners],
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
