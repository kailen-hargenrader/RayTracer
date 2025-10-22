import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# specify cube data

data = {
    "MESH": {
        "cube": {
            1: {"translation": [0, 0, 0], "rotation": [0, 0, 0], "scale": [1.0, 0.5, 1.0]},
            2: {"translation": [2, 0, 0], "rotation": [0, 0.5, 0], "scale": [1.0, 0.5, 1.5]},
            3: {"translation": [0, 2, 0], "rotation": [0.3, 0.6, 0], "scale": [.2, 1.5, 1.5]},
        }
    }
}

# Load cubes into Blender

def load_cubes_from_json():
	cubes = data.get("MESH", {}).get("cube", {})
	for _id, params in cubes.items():
		trans = params.get("translation", [0, 0, 0])
		rot = params.get("rotation", [0, 0, 0])
		scale = params.get("scale", [1.0, 1.0, 1.0])
		bpy.ops.mesh.primitive_cube_add()
		obj = bpy.context.active_object
		obj.location = trans
		obj.rotation_euler = rot
		obj.scale = tuple(scale)

if __name__ == "__main__":
	load_cubes_from_json()
