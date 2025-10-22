import json
import math
import os
import sys

try:
	bpy = __import__("bpy")
	from mathutils import Vector
except Exception:
	raise SystemExit("Run inside Blender")

# Optional: used for plane mesh creation
try:
	bmesh = __import__("bmesh")
except Exception:
	bmesh = None


def _as_vec3(v, default=(0.0, 0.0, 0.0)):
	if v is None:
		return default
	if isinstance(v, (list, tuple)) and len(v) >= 3:
		return (float(v[0]), float(v[1]), float(v[2]))
	if isinstance(v, dict):
		x = float(v.get("x", default[0]))
		y = float(v.get("y", default[1]))
		z = float(v.get("z", default[2]))
		return (x, y, z)
	return default


def _as_vec2(v, default=(1920, 1080)):
	if v is None:
		return default
	if isinstance(v, (list, tuple)):
		if len(v) >= 2:
			return (int(v[0]), int(v[1]))
		return default
	if isinstance(v, dict):
		x_val = v.get("x")
		y_val = v.get("y")
		if x_val is None or y_val is None:
			if default is None:
				return None
			return (int(default[0]), int(default[1]))
		return (int(x_val), int(y_val))
	return default


def _as_euler(rot):
	# Accept either {roll,pitch,yaw} or [r,p,y]
	if isinstance(rot, dict):
		return (
			float(rot.get("roll", 0.0)),
			float(rot.get("pitch", 0.0)),
			float(rot.get("yaw", 0.0)),
		)
	if isinstance(rot, (list, tuple)) and len(rot) >= 3:
		return (float(rot[0]), float(rot[1]), float(rot[2]))
	return (0.0, 0.0, 0.0)


def _get_first(container, keys, default=None):
	for k in keys:
		if k in container:
			return container[k]
	return default if default is not None else {}


def _load_cubes(mesh_section):
	cubes = _get_first(mesh_section, ["Cube", "cube"], {})
	for _id, p in (cubes or {}).items():
		trans = _as_vec3(p.get("translation"))
		rot = _as_euler(p.get("rotation"))
		scale_v = _as_vec3(p.get("scale"), (1.0, 1.0, 1.0))
		bpy.ops.mesh.primitive_cube_add()
		obj = bpy.context.active_object
		obj.location = trans
		obj.rotation_euler = rot
		obj.scale = scale_v


def _load_cylinders(mesh_section):
	cyls = _get_first(mesh_section, ["Cylinder", "cylinder"], {})
	for _id, p in (cyls or {}).items():
		trans = _as_vec3(p.get("translation"))
		rot = _as_euler(p.get("rotation"))
		scale_v = _as_vec3(p.get("scale"), (1.0, 1.0, 1.0))
		bpy.ops.mesh.primitive_cylinder_add()
		obj = bpy.context.active_object
		obj.location = trans
		obj.rotation_euler = rot
		obj.scale = scale_v


def _load_spheres(mesh_section):
	spheres = _get_first(mesh_section, ["Sphere", "sphere"], {})
	for _id, p in (spheres or {}).items():
		loc = _as_vec3(p.get("location"))
		scale_v = _as_vec3(p.get("scale"), (1.0, 1.0, 1.0))
		bpy.ops.mesh.primitive_uv_sphere_add()
		obj = bpy.context.active_object
		obj.location = loc
		obj.scale = scale_v


def _load_planes(mesh_section):
	planes = _get_first(mesh_section, ["Plane", "plane"], {})
	for _id, p in (planes or {}).items():
		corners = p.get("corners", [])
		if len(corners) != 4:
			continue
		coords = [_as_vec3(c) for c in corners]
		if bmesh is None:
			# Fallback: add a default plane
			bpy.ops.mesh.primitive_plane_add()
			obj = bpy.context.active_object
			# Without bmesh, we won't reshape to the exact quad
			continue
		mesh = bpy.data.meshes.new("Plane")
		bm = bmesh.new()
		verts = [bm.verts.new(co) for co in coords]
		bm.faces.new(verts)
		bm.to_mesh(mesh)
		bm.free()
		obj = bpy.data.objects.new("Plane", mesh)
		bpy.context.collection.objects.link(obj)


def _load_point_lights(light_section):
	pts = _get_first(light_section, ["POINT", "point"], {})
	for _id, p in (pts or {}).items():
		loc = _as_vec3(p.get("location"))
		ri = float(p.get("radiant_intensity", 1.0))
		bpy.ops.object.light_add(type='POINT')
		light = bpy.context.active_object
		light.location = loc
		light.data.energy = 4.0 * math.pi * ri


def _load_persp_cameras(cam_section):
	per = _get_first(cam_section, ["PERSP", "persp"], {})
	for _id, p in (per or {}).items():
		loc = _as_vec3(p.get("location"))
		dirv = _as_vec3(p.get("direction"), (0.0, 0.0, -1.0))
		fl = float(p.get("focal_length", 50.0))
		sw = float(p.get("sensor_width", 36.0))
		sh = float(p.get("sensor_height", 24.0))
		res = _as_vec2(p.get("film_resolution"), None)
		bpy.ops.object.camera_add()
		cam = bpy.context.active_object
		cam.location = loc
		# Orient camera so its local -Z points along the desired world direction
		dir_vec = Vector(dirv)
		if dir_vec.length == 0:
			dir_vec = Vector((0.0, 0.0, -1.0))
		q = dir_vec.normalized().to_track_quat('-Z', 'Y')
		cam.rotation_euler = q.to_euler('XYZ')
		cam.data.type = 'PERSP'
		cam.data.lens = fl
		cam.data.sensor_width = sw
		cam.data.sensor_height = sh
		if res:
			scene = bpy.context.scene
			scene.render.resolution_x = int(res[0])
			scene.render.resolution_y = int(res[1])


def load_from_json(path):
	with open(path, "r") as f:
		data = json.load(f)

	mesh = _get_first(data, ["MESH", "mesh"], {})
	light = _get_first(data, ["LIGHT", "light"], {})
	cam = _get_first(data, ["CAMERA", "camera"], {})

	if mesh:
		_load_cubes(mesh)
		_load_spheres(mesh)
		_load_planes(mesh)
		_load_cylinders(mesh)
	if light:
		_load_point_lights(light)
	if cam:
		_load_persp_cameras(cam)


def _resolve_default_path():
	# Prefer CLI arg after '--'; else try common defaults
	argv = sys.argv
	if "--" in argv:
		idx = argv.index("--")
		if idx + 1 < len(argv):
			return argv[idx + 1]
	# Try repo ASCII path then local
	path = os.path.join(os.path.expanduser("~"), "RayTracer", "s2899221", "ASCII", "cube_camera.json")
	if os.path.exists(path):
		return path
	else:
		raise ValueError(f"File not found: {path}")


if __name__ == "__main__":
	json_path = _resolve_default_path()
	print(f"Loading scene from JSON: {json_path}")
	load_from_json(json_path)
	print("Done.")


