import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# Write sample JSON with multiple spheres

data = {
    "MESH": {
        "sphere": {
            1: {"location": [0, 0, 1], "radius": 1.0},
            2: {"location": [2, 0, 1], "radius": 0.5},
            3: {"location": [0, 2, 1], "radius": 1.5},
        }
    }
}
data["MESH"]["sphere"] = {str(k): v for k, v in data["MESH"]["sphere"].items()}

# Load spheres into Blender

def load_spheres_from_json(path="spheres.json"):
	spheres = data.get("MESH", {}).get("sphere", {})
	for _id, params in spheres.items():
		loc = params.get("location", [0, 0, 0])
		radius = float(params.get("radius", 1.0))
		bpy.ops.mesh.primitive_uv_sphere_add()
		obj = bpy.context.active_object
		obj.location = loc
		obj.scale = (radius, radius, radius)

if __name__ == "__main__":
	load_spheres_from_json()
