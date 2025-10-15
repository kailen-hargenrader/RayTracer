import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# specify cylinder data

data = {
	"MESH": {
		"cylinder": {
			1: {"translation": [0, 0, 0], "rotation": [0, 0, 0], "radius": 1.0, "length": 2.0},
			2: {"translation": [2, 0, 0], "rotation": [0, 0.5, 0], "radius": 0.5, "length": 1.0},
			3: {"translation": [0, 2, 0], "rotation": [0.3, 0.6, 0], "radius": 1.5, "length": 3.0},
		}
	}
}

# Load cylinders into Blender

def load_cylinders_from_json():
	cylinders = data.get("MESH", {}).get("cylinder", {})
	for _id, params in cylinders.items():
		trans = params.get("translation", [0, 0, 0])
		rot = params.get("rotation", [0, 0, 0])
		radius = float(params.get("radius", 1.0))
		length = float(params.get("length", 2.0))
		bpy.ops.mesh.primitive_cylinder_add()
		obj = bpy.context.active_object
		obj.location = trans
		obj.rotation_euler = rot
		# Blender's default cylinder has radius=1 and depth=2 along Z
		obj.scale = (radius, radius, length / 2.0)

if __name__ == "__main__":
	load_cylinders_from_json()


