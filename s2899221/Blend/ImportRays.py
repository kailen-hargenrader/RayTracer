import bpy
from mathutils import Vector

# ===== Configuration =====
# Change this to the path of your generated rays file (out.txt)
FILE_PATH = r"C:\\Users\\Kailen\\RayTracer\\out.txt"

RAY_LENGTH = 1.0           # How long to draw each ray
CYLINDER_RADIUS = 0.02     # Radius of the cylinder representing the ray
SPHERE_RADIUS = 0.04       # Radius of the sphere at the ray tip

CREATE_COLLECTION = True
COLLECTION_NAME = "RaysFromFile"


def parse_rays(filepath):
    rays = []  # list of (origin Vector, direction Vector)
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                ox, oy, oz, dx, dy, dz = map(float, parts[:6])
            except ValueError:
                continue
            origin = Vector((ox, oy, oz))
            direction = Vector((dx, dy, dz))
            if direction.length == 0.0:
                continue
            direction.normalize()
            rays.append((origin, direction))
    return rays


def ensure_collection(name):
    coll = bpy.data.collections.get(name)
    if coll is None:
        coll = bpy.data.collections.new(name)
        bpy.context.scene.collection.children.link(coll)
    return coll


def link_to_collection(obj, collection):
    if collection is None:
        return
    for c in obj.users_collection:
        c.objects.unlink(obj)
    collection.objects.link(obj)


def create_ray_objects(origin, direction, length, coll=None, index=0):
    end = origin + direction * length
    mid = origin + (direction * (length * 0.5))

    # Create cylinder aligned with +Z by default
    bpy.ops.mesh.primitive_cylinder_add(radius=CYLINDER_RADIUS, depth=length, location=mid)
    cyl = bpy.context.active_object
    cyl.name = f"Ray_{index:04d}"

    # Rotate cylinder so local +Z aligns with the ray direction
    rot_quat = Vector((0.0, 0.0, 1.0)).rotation_difference(direction)
    cyl.rotation_euler = rot_quat.to_euler()

    # Move to collection if needed
    link_to_collection(cyl, coll)

    # Create sphere at the end to indicate direction
    bpy.ops.mesh.primitive_uv_sphere_add(radius=SPHERE_RADIUS, location=end)
    sph = bpy.context.active_object
    sph.name = f"RayTip_{index:04d}"
    link_to_collection(sph, coll)

    # Parent the sphere to the cylinder while preserving world transform
    sph.parent = cyl
    sph.matrix_parent_inverse = cyl.matrix_world.inverted()


def main():
    rays = parse_rays(FILE_PATH)
    if not rays:
        print("No rays parsed from file. Check FILE_PATH and file contents.")
        return

    coll = ensure_collection(COLLECTION_NAME) if CREATE_COLLECTION else None

    # Ensure we're in Object mode
    if bpy.ops.object.mode_set.poll():
        try:
            bpy.ops.object.mode_set(mode='OBJECT')
        except Exception:
            pass

    for i, (origin, direction) in enumerate(rays):
        create_ray_objects(origin, direction, RAY_LENGTH, coll=coll, index=i)

    print(f"Created {len(rays)} rays from '{FILE_PATH}'.")


if __name__ == "__main__":
    main()


