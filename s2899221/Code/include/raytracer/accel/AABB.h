#pragma once

#include "raytracer/math/Vec3.h"
#include "raytracer/utils/Ray.h"
#include <algorithm>

namespace rt {

struct AABB {
    Vec3 min { 1e30f, 1e30f, 1e30f };
    Vec3 max { -1e30f, -1e30f, -1e30f };

    AABB() = default;
    AABB(const Vec3& mn, const Vec3& mx) : min(mn), max(mx) {}

    static AABB fromPoints(const Vec3& a, const Vec3& b) {
        return AABB(Vec3::min(a,b), Vec3::max(a,b));
    }
    void expand(const Vec3& p) {
        min = Vec3::min(min, p);
        max = Vec3::max(max, p);
    }
    void expand(const AABB& b) {
        expand(b.min);
        expand(b.max);
    }
    bool isValid() const { return min.x <= max.x && min.y <= max.y && min.z <= max.z; }

    bool intersect(const Ray& r, float tMin, float tMax) const {
        for (int axis = 0; axis < 3; ++axis) {
            float ro = r.origin[axis];
            float rd = r.direction[axis];
            float inv = (std::abs(rd) > 1e-8f) ? 1.0f / rd : 1e30f;
            float t0 = (min[axis] - ro) * inv;
            float t1 = (max[axis] - ro) * inv;
            if (t0 > t1) std::swap(t0, t1);
            tMin = t0 > tMin ? t0 : tMin;
            tMax = t1 < tMax ? t1 : tMax;
            if (tMax < tMin) return false;
        }
        return true;
    }
};

AABB transformAABB(const AABB& box, const class Mat4& m);

} // namespace rt


