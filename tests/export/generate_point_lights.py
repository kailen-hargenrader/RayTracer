import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# Write sample JSON with multiple point lights

def write_point_lights_json(path="point_lights.json"):
	data = {
		"LIGHT": {
			"point": {
				1: {"location": [0, 0, 3], "radiant_intensity": 10.0},
				2: {"location": [3, 0, 3], "radiant_intensity": 5.0},
				3: {"location": [0, 3, 3], "radiant_intensity": 20.0},
			}
		}
	}
	data["LIGHT"]["point"] = {str(k): v for k, v in data["LIGHT"]["point"].items()}
	with open(path, "w") as f:
		json.dump(data, f, indent=4)

# Load point lights into Blender

def load_point_lights_from_json(path="point_lights.json"):
	with open(path, "r") as f:
		data = json.load(f)
	points = data.get("LIGHT", {}).get("point", {})
	for _id, params in points.items():
		loc = params.get("location", [0, 0, 0])
		ri = float(params.get("radiant_intensity", 1.0))
		bpy.ops.object.light_add(type='POINT')
		light = bpy.context.active_object
		light.location = loc
		# Convert radiant intensity back to approximate energy: energy = 4*pi*I
		light.data.energy = 4 * 3.141592653589793 * ri

if __name__ == "__main__":
	write_point_lights_json()
	load_point_lights_from_json()
