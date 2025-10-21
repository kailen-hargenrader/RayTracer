import bpy
import json
import random
from math import sqrt
from mathutils import Vector


# ============ Configuration ============
# Path to scene JSON (meshes/cameras/lights). Update as needed.
JSON_PATH = r"C:\\Users\\Kailen\\RayTracer\\s2899221\\ASCII\\lab2.json"

# Number of random rays to sample from the first perspective camera
NUM_RAYS = 200

# Visual parameters
RAY_CYL_RADIUS = 0.015
REFLECT_CYL_RADIUS = 0.012
MAX_RAY_LENGTH = 5.0
REFLECT_LENGTH = 2.0

# Collection names
RAYS_COLL_NAME = "RayIntersections"
REFL_COLL_NAME = "RayReflections"


# ============ Helpers ============
def ensure_collection(name: str):
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


def add_cylinder_between(p0: Vector, p1: Vector, radius: float, name: str, collection=None):
    direction = p1 - p0
    length = direction.length
    if length <= 1e-9:
        return None
    mid = p0 + direction * 0.5
    direction.normalize()

    bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=length, location=mid)
    cyl = bpy.context.active_object
    cyl.name = name
    rot_quat = Vector((0.0, 0.0, 1.0)).rotation_difference(direction)
    cyl.rotation_euler = rot_quat.to_euler()
    link_to_collection(cyl, collection)
    return cyl


def build_camera_basis(location: Vector, direction: Vector):
    forward = direction.normalized()
    world_up = Vector((0.0, 0.0, 1.0))
    right = forward.cross(world_up).normalized()
    if right.length < 1e-6:
        world_up = Vector((0.0, 1.0, 0.0))
        right = forward.cross(world_up).normalized()
    up = right.cross(forward)
    return right, up, forward


def pixel_to_ray(px: float, py: float, cam):
    # cam: dict with fields 'position' (Vector), 'right','up','forward' (Vector),
    #      'focal','sensor_w','sensor_h','resx','resy'
    nx = px / float(cam['resx'])
    ny = py / float(cam['resy'])
    x_sensor = cam['sensor_w'] * (nx - 0.5)
    y_sensor = cam['sensor_h'] * (0.5 - ny)
    dir_world = (
        cam['right'] * x_sensor +
        cam['up'] * y_sensor +
        cam['forward'] * cam['focal']
    )
    return cam['position'], dir_world.normalized()


# ============ Intersections ============
def rpy_to_matrix(roll: float, pitch: float, yaw: float):
    # Construct R = Rz(yaw) * Ry(pitch) * Rx(roll)
    import math
    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)
    Rz = ((cy, -sy, 0.0), (sy, cy, 0.0), (0.0, 0.0, 1.0))
    Ry = ((cp, 0.0, sp), (0.0, 1.0, 0.0), (-sp, 0.0, cp))
    Rx = ((1.0, 0.0, 0.0), (0.0, cr, -sr), (0.0, sr, cr))
    # A = Ry * Rx
    A = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            A[i][j] = Ry[i][0]*Rx[0][j] + Ry[i][1]*Rx[1][j] + Ry[i][2]*Rx[2][j]
    # R = Rz * A
    R = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            R[i][j] = Rz[i][0]*A[0][j] + Rz[i][1]*A[1][j] + Rz[i][2]*A[2][j]
    return R


def mul_mat_vec(M, v: Vector) -> Vector:
    return Vector((
        M[0][0]*v.x + M[0][1]*v.y + M[0][2]*v.z,
        M[1][0]*v.x + M[1][1]*v.y + M[1][2]*v.z,
        M[2][0]*v.x + M[2][1]*v.y + M[2][2]*v.z,
    ))


def transpose3(M):
    return (
        (M[0][0], M[1][0], M[2][0]),
        (M[0][1], M[1][1], M[2][1]),
        (M[0][2], M[1][2], M[2][2]),
    )


def reflect(dir_world: Vector, n_world: Vector) -> Vector:
    n_unit = n_world.normalized()
    return (dir_world - 2.0 * dir_world.dot(n_unit) * n_unit).normalized()


def intersect_cube(origin: Vector, direction: Vector, tr: Vector, rpy, scale: float):
    # Build transforms
    R = rpy_to_matrix(rpy[0], rpy[1], rpy[2])
    RT = transpose3(R)
    # to object frame
    to = origin - tr
    ro_r = mul_mat_vec(RT, to)
    rd_r = mul_mat_vec(RT, direction)
    invS = 1.0/scale if scale != 0.0 else 0.0
    ro = Vector((ro_r.x*invS, ro_r.y*invS, ro_r.z*invS))
    rd = Vector((rd_r.x*invS, rd_r.y*invS, rd_r.z*invS))

    eps = 1e-12
    tmin = -1e30
    tmax = 1e30
    n_enter = Vector((0.0,0.0,0.0))
    n_exit  = Vector((0.0,0.0,0.0))
    for axis in range(3):
        o = ro[axis]
        d = rd[axis]
        low, high = -1.0, 1.0
        if abs(d) < eps:
            if o < low or o > high:
                return (False, None)
            continue
        t1 = (low - o)/d
        t2 = (high - o)/d
        n1 = Vector((0.0,0.0,0.0)); n1[axis] = -1.0
        n2 = Vector((0.0,0.0,0.0)); n2[axis] =  1.0
        if t1 > t2:
            t1, t2 = t2, t1
            n1, n2 = n2, n1
        if t1 > tmin:
            tmin = t1; n_enter = n1
        if t2 < tmax:
            tmax = t2; n_exit = n2
        if tmin > tmax:
            return (False, None)
    if tmax < 0.0:
        return (False, None)
    t_hit = tmin if tmin >= 0.0 else tmax
    n_local = n_enter if tmin >= 0.0 else n_exit

    p_local = ro + rd * t_hit
    # world normal and point
    # n_world = R * n_local
    n_world = mul_mat_vec(R, n_local).normalized()
    p_world = mul_mat_vec(R, Vector((p_local.x*scale, p_local.y*scale, p_local.z*scale))) + tr
    dist = (p_world - origin).length
    refl = reflect(direction, n_world)
    return (True, {
        'distance': dist,
        'point': p_world,
        'normal': n_world,
        'reflection': refl,
    })


def point_in_quad(p: Vector, corners: list, n: Vector) -> bool:
    # All edges must have the same sign of dot(n, cross(edge, p-vi))
    signs = []
    for i in range(4):
        vi = corners[i]
        vj = corners[(i+1) % 4]
        edge = vj - vi
        cp = edge.cross(p - vi)
        signs.append(n.dot(cp) >= -1e-9)
    return all(signs) or not any(signs)


def intersect_plane(origin: Vector, direction: Vector, corners: list):
    if len(corners) < 3:
        return (False, None)
    v0, v1, v3 = corners[0], corners[1], corners[3]
    n = (v1 - v0).cross(v3 - v0)
    if n.length < 1e-12:
        return (False, None)
    n = n.normalized()
    denom = n.dot(direction)
    if abs(denom) < 1e-12:
        return (False, None)
    t = n.dot(v0 - origin) / denom
    if t < 0.0:
        return (False, None)
    p = origin + direction * t
    if not point_in_quad(p, corners, n):
        return (False, None)
    refl = reflect(direction, n)
    return (True, {
        'distance': (p - origin).length,
        'point': p,
        'normal': n,
        'reflection': refl,
    })


# ============ JSON Parsing ============
def parse_scene(json_path: str):
    with open(json_path, 'r') as f:
        data = json.load(f)

    # Camera (first perspective camera)
    cam_block = data.get('CAMERA', {}).get('PERSP', {})
    if not cam_block:
        raise RuntimeError('No CAMERA.PERSP found in JSON')
    first_cam_id = next(iter(cam_block))
    c = cam_block[first_cam_id]
    cam_location = Vector((c['location']['x'], c['location']['y'], c['location']['z']))
    cam_direction = Vector((c['direction']['x'], c['direction']['y'], c['direction']['z']))
    right, up, forward = build_camera_basis(cam_location, cam_direction)
    camera = {
        'position': cam_location,
        'right': right,
        'up': up,
        'forward': forward,
        'focal': float(c['focal_length']),
        'sensor_w': float(c['sensor_width']),
        'sensor_h': float(c['sensor_height']),
        'resx': int(c['film_resolution']['x']),
        'resy': int(c['film_resolution']['y']),
    }

    # Meshes
    meshes = []
    mblock = data.get('MESH', {})
    # Cubes
    for _id, cube in mblock.get('Cube', {}).items():
        tr = cube.get('translation', {'x':0,'y':0,'z':0})
        rot = cube.get('rotation', {'roll':0,'pitch':0,'yaw':0})
        sc = float(cube.get('scale', 1.0))
        meshes.append({
            'type': 'Cube',
            'translation': Vector((float(tr['x']), float(tr['y']), float(tr['z']))),
            'rotation': (float(rot['roll']), float(rot['pitch']), float(rot['yaw'])),
            'scale': sc,
        })
    # Planes
    for _id, plane in mblock.get('Plane', {}).items():
        crs = []
        for pt in plane.get('corners', []):
            crs.append(Vector((float(pt['x']), float(pt['y']), float(pt['z']))))
        if len(crs) == 4:
            meshes.append({ 'type': 'Plane', 'corners': crs })

    return camera, meshes


# ============ Main ============
def main(json_path: str = JSON_PATH, num_rays: int = NUM_RAYS):
    camera, meshes = parse_scene(json_path)

    rays_coll = ensure_collection(RAYS_COLL_NAME)
    refl_coll = ensure_collection(REFL_COLL_NAME)

    # Ensure Object mode
    if bpy.ops.object.mode_set.poll():
        try:
            bpy.ops.object.mode_set(mode='OBJECT')
        except Exception:
            pass

    # Sample random pixels
    for i in range(num_rays):
        px = random.uniform(0.0, float(camera['resx']))
        py = random.uniform(0.0, float(camera['resy']))
        origin, direction = pixel_to_ray(px, py, camera)

        # Intersect with all meshes
        best_hit = None
        best_dist = 1e30
        for m in meshes:
            if m['type'] == 'Cube':
                ok, info = intersect_cube(origin, direction, m['translation'], m['rotation'], m['scale'])
            elif m['type'] == 'Plane':
                ok, info = intersect_plane(origin, direction, m['corners'])
            else:
                ok, info = (False, None)
            if ok and info['distance'] < best_dist:
                best_dist = info['distance']
                best_hit = info

        # Draw primary ray
        if best_hit is not None:
            end = best_hit['point']
        else:
            end = origin + direction * MAX_RAY_LENGTH
        add_cylinder_between(origin, end, RAY_CYL_RADIUS, name=f"Ray_{i:04d}", collection=rays_coll)

        # Draw reflection ray if hit
        if best_hit is not None:
            refl_end = best_hit['point'] + best_hit['reflection'] * REFLECT_LENGTH
            add_cylinder_between(best_hit['point'], refl_end, REFLECT_CYL_RADIUS, name=f"Refl_{i:04d}", collection=refl_coll)

    print(f"Generated {num_rays} rays and reflections from '{json_path}'.")


if __name__ == "__main__":
    main()


