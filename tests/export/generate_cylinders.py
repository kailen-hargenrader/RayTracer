import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# specify cylinder data

data = {
	"MESH": {
		"cylinder": {
			1: {"translation": [0, 0, 0], "rotation": [0, 0, 0], "scale": [.5, 1.0, 1.5]},
			2: {"translation": [2, 0, 0], "rotation": [0, 0.5, 0], "scale": [0.2, 0.4, 0.6]},
			3: {"translation": [0, 2, 0], "rotation": [0.3, 0.6, 0], "scale": [1.5, 1.5, 1.5]},
		}
	}
}

# Load cylinders into Blender

def load_cylinders_from_json():
	cylinders = data.get("MESH", {}).get("cylinder", {})
	for _id, params in cylinders.items():
		trans = params.get("translation", [0, 0, 0])
		rot = params.get("rotation", [0, 0, 0])
		scale = params.get("scale", [1.0, 1.0, 1.0])
		bpy.ops.mesh.primitive_cylinder_add()
		obj = bpy.context.active_object
		obj.location = trans
		obj.rotation_euler = rot
		obj.scale = tuple(scale)

if __name__ == "__main__":
	load_cylinders_from_json()


