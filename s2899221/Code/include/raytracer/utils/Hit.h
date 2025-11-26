#pragma once

#include <cfloat>
#include <cstdint>
#include "raytracer/math/Vec2.h"
#include "raytracer/math/Vec3.h"
#include "raytracer/utils/Ray.h"

namespace rt {

struct Hit {
    Ray ray;
    bool didHit {false};
    float t {FLT_MAX};
    Vec3 position;     // world-space
    Vec3 normal;       // world-space, unit-length, facing against ray if outside
    Vec2 uv;           // optional surface UV if available

    // Simple shading fields (extend later as needed)
    Vec3 albedo {1.0f, 1.0f, 1.0f};
    bool frontFace {true};

    void setFaceNormal(const Vec3& outwardNormal) {
        frontFace = Vec3::dot(ray.direction, outwardNormal) < 0.0f;
        normal = frontFace ? outwardNormal : outwardNormal * -1.0f;
    }
};

} // namespace rt


