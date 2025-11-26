#pragma once

#include "raytracer/math/Vec3.h"

namespace rt {

struct Ray {
    Vec3 origin;
    Vec3 direction; // Should be normalized

    Ray() = default;
    Ray(const Vec3& o, const Vec3& d) : origin(o), direction(d) {}

    Vec3 at(float t) const { return origin + direction * t; }
};

} // namespace rt


