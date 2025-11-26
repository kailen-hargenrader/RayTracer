#include "raytracer/accel/AABB.h"
#include "raytracer/math/Mat4.h"

namespace rt {

static inline Vec3 transformPointSafe(const Mat4& m, const Vec3& p) {
    return m.transformPoint(p);
}

AABB transformAABB(const AABB& box, const Mat4& m) {
    // Transform 8 corners and recompute bounds
    Vec3 corners[8] = {
        {box.min.x, box.min.y, box.min.z},
        {box.max.x, box.min.y, box.min.z},
        {box.min.x, box.max.y, box.min.z},
        {box.min.x, box.min.y, box.max.z},
        {box.max.x, box.max.y, box.min.z},
        {box.max.x, box.min.y, box.max.z},
        {box.min.x, box.max.y, box.max.z},
        {box.max.x, box.max.y, box.max.z},
    };
    AABB out;
    for (int i = 0; i < 8; ++i) {
        out.expand(transformPointSafe(m, corners[i]));
    }
    return out;
}

} // namespace rt


