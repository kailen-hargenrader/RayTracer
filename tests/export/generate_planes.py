import json
try:
	import bpy
except ImportError:
	raise SystemExit("Run inside Blender")


# Define two example quads
p1 = {
    "corners": [[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]]
}
p2 = {
    "corners": [[0, 0, 0], [2, 0, 0], [2, 0, 2], [0, 0, 2]]
}
data = {"MESH": {"plane": {1: p1, 2: p2}}}
data["MESH"]["plane"] = {str(k): v for k, v in data["MESH"]["plane"].items()}

# Load planes into Blender using bmesh to build quad from corners

def load_planes_from_json(path="planes.json"):
	import bmesh
	planes = data.get("MESH", {}).get("plane", {})
	for _id, params in planes.items():
		corners = params.get("corners", [])
		if len(corners) != 4:
			continue
		mesh = bpy.data.meshes.new("Plane")
		bm = bmesh.new()
		verts = [bm.verts.new(tuple(c)) for c in corners]
		bm.faces.new(verts)
		bm.to_mesh(mesh)
		bm.free()
		obj = bpy.data.objects.new("Plane", mesh)
		bpy.context.collection.objects.link(obj)

if __name__ == "__main__":
	load_planes_from_json()
