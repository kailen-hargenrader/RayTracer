import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# Optional: used for plane mesh creation
try:
	import bmesh
except Exception:
	bmesh = None


def vec3(x, y, z):
	return {"x": float(x), "y": float(y), "z": float(z)}


def vec2(x, y):
	return {"x": int(x), "y": int(y)}


def dict3_to_tuple(d):
	return (float(d.get("x", 0.0)), float(d.get("y", 0.0)), float(d.get("z", 0.0)))


data = {
	"MESH": {
		"Cube": {
			"1": {"translation": vec3(-2.0, 1.5, 0.5), "rotation": {"roll": 0.25, "pitch": 0.5, "yaw": 0.75}, "scale": 1.0},
			"2": {"translation": vec3(3.0, -1.0, 1.0), "rotation": {"roll": 1.1, "pitch": 0.2, "yaw": 0.0}, "scale": 0.65},
			"3": {"translation": vec3(0.0, 3.0, 2.0), "rotation": {"roll": 0.0, "pitch": 0.9, "yaw": 1.2}, "scale": 1.8},
			"4": {"translation": vec3(-3.5, -2.0, 0.2), "rotation": {"roll": 0.4, "pitch": 0.8, "yaw": 0.1}, "scale": 1.25},
		},
		"Sphere": {
			"1": {"location": vec3(1.0, 1.0, 1.0), "radius": 1.0},
			"2": {"location": vec3(-1.5, 2.5, 0.5), "radius": 0.6},
			"3": {"location": vec3(2.0, -2.5, 1.5), "radius": 1.4},
			"4": {"location": vec3(-2.0, -1.0, 2.0), "radius": 0.9},
		},
		"Plane": {
			"1": {"corners": [vec3(-2, -2, 0), vec3(2, -2, 0), vec3(2, 2, 0), vec3(-2, 2, 0)]},
			"2": {"corners": [vec3(0, 0, 1), vec3(3, 0, 1), vec3(3, 2, 1), vec3(0, 2, 1)]},
			"3": {"corners": [vec3(-1, 3, 0.5), vec3(1, 3, 0.5), vec3(1, 5, 0.5), vec3(-1, 5, 0.5)]},
		},
	},
	"LIGHT": {
		"POINT": {
			"1": {"location": vec3(0.0, 0.0, 4.0), "radiant_intensity": 12.0},
			"2": {"location": vec3(4.0, -3.0, 3.0), "radiant_intensity": 6.5},
			"3": {"location": vec3(-4.0, 2.0, 5.0), "radiant_intensity": 18.0},
		},
	},
	"CAMERA": {
		"PERSP": {
			"1": {
				"location": vec3(-6.0, -6.0, 3.0),
				"direction": vec3(0.6, 0.7, -0.2),
				"focal_length": 50.0,
				"sensor_width": 36.0,
				"sensor_height": 24.0,
				"film_resolution": vec2(1920, 1080),
			},
			"2": {
				"location": vec3(8.0, -5.0, 4.0),
				"direction": vec3(-0.8, 0.5, -0.1),
				"focal_length": 35.0,
				"sensor_width": 36.0,
				"sensor_height": 24.0,
				"film_resolution": vec2(1280, 720),
			},
		},
	},
}





# Loader utilities


def load_from_json(path="random.json"):

	# Load cubes
	for _id, p in data.get("MESH", {}).get("Cube", {}).items():
		trans = dict3_to_tuple(p.get("translation", {}))
		rot = p.get("rotation", {"roll": 0.0, "pitch": 0.0, "yaw": 0.0})
		scale = float(p.get("scale", 1.0))
		bpy.ops.mesh.primitive_cube_add()
		obj = bpy.context.active_object
		obj.location = trans
		obj.rotation_euler = (float(rot.get("roll", 0.0)), float(rot.get("pitch", 0.0)), float(rot.get("yaw", 0.0)))
		obj.scale = (scale, scale, scale)

	# Load spheres
	for _id, p in data.get("MESH", {}).get("Sphere", {}).items():
		loc = dict3_to_tuple(p.get("location", {}))
		radius = float(p.get("radius", 1.0))
		bpy.ops.mesh.primitive_uv_sphere_add()
		obj = bpy.context.active_object
		obj.location = loc
		obj.scale = (radius, radius, radius)

	# Load planes
	planes = data.get("MESH", {}).get("Plane", {})
	for _id, p in planes.items():
		corners = p.get("corners", [])
		if len(corners) != 4:
			continue
		c = [dict3_to_tuple(corner) for corner in corners]
		if bmesh is None:
			# Fallback: add a default plane and move vertices in edit mode
			bpy.ops.mesh.primitive_plane_add()
			obj = bpy.context.active_object
			# Can't easily set arbitrary quad without bmesh; skip precise placement
			continue
		mesh = bpy.data.meshes.new("Plane")
		bm = bmesh.new()
		verts = [bm.verts.new(coord) for coord in c]
		bm.faces.new(verts)
		bm.to_mesh(mesh)
		bm.free()
		obj = bpy.data.objects.new("Plane", mesh)
		bpy.context.collection.objects.link(obj)

	# Load point lights
	for _id, p in data.get("LIGHT", {}).get("POINT", {}).items():
		loc = dict3_to_tuple(p.get("location", {}))
		ri = float(p.get("radiant_intensity", 1.0))
		bpy.ops.object.light_add(type='POINT')
		light = bpy.context.active_object
		light.location = loc
		light.data.energy = 4 * 3.141592653589793 * ri

	# Load perspective cameras
	for _id, p in data.get("CAMERA", {}).get("PERSP", {}).items():
		loc = dict3_to_tuple(p.get("location", {}))
		dirv = dict3_to_tuple(p.get("direction", {}))
		fl = float(p.get("focal_length", 50.0))
		sw = float(p.get("sensor_width", 36.0))
		sh = float(p.get("sensor_height", 24.0))
		bpy.ops.object.camera_add()
		cam = bpy.context.active_object
		cam.location = loc
		# Simple orientation: compute yaw/pitch from direction and set rotation_euler
		import math
		x, y, z = dirv
		# Avoid zero-length
		l = math.sqrt(x*x + y*y + z*z) or 1.0
		x /= l; y /= l; z /= l
		# Yaw around Z, pitch around X to roughly align -Z to dir
		yaw = math.atan2(y, x)
		pitch = -math.asin(z)
		cam.rotation_euler = (pitch, 0.0, yaw)
		cam.data.type = 'PERSP'
		cam.data.lens = fl
		cam.data.sensor_width = sw
		cam.data.sensor_height = sh

if __name__ == "__main__":
	load_from_json()
