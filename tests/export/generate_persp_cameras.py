import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")

# Write sample JSON with multiple perspective cameras

data = {
    "CAMERA": {
        "PERSP": {
            1: {
                "location": [0, -5, 2],
                "direction": [0, 1, -0.3],
                "focal_length": 50.0,
                "sensor_width": 36.0,
                "sensor_height": 24.0,
                "film_resolution": [1920, 1080],
            },
            2: {
                "location": [5, -5, 3],
                "direction": [-1, 1, -0.2],
                "focal_length": 35.0,
                "sensor_width": 36.0,
                "sensor_height": 24.0,
                "film_resolution": [1280, 720],
            },
        }
    }
}
data["CAMERA"]["PERSP"] = {str(k): v for k, v in data["CAMERA"]["PERSP"].items()}

# Load perspective cameras into Blender

def load_persp_cameras_from_json(path="persp_cameras.json"):
	import mathutils
	persp = data.get("CAMERA", {}).get("PERSP", {})
	for _id, params in persp.items():
		loc = params.get("location", [0, 0, 0])
		dirv = params.get("direction", [0, 0, -1])
		fl = float(params.get("focal_length", 50.0))
		sw = float(params.get("sensor_width", 36.0))
		sh = float(params.get("sensor_height", 24.0))
		bpy.ops.object.camera_add()
		cam = bpy.context.active_object
		cam.location = loc
		# Orient camera so its -Z aligns with given direction
		dir_vec = mathutils.Vector(dirv).normalized()
		# Compute rotation: track -Z to dir, up Y
		track = -dir_vec
		up = mathutils.Vector((0, 1, 0))
		rot = track.to_track_quat('Z', 'Y').to_euler()
		cam.rotation_euler = rot
		cam.data.type = 'PERSP'
		cam.data.lens = fl
		cam.data.sensor_width = sw
		cam.data.sensor_height = sh

if __name__ == "__main__":
	load_persp_cameras_from_json()
