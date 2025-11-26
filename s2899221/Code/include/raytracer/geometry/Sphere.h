#pragma once

#include "raytracer/geometry/Mesh.h"

namespace rt {

// Canonical unit sphere at origin with radius 0.5 in object space
class Sphere : public Mesh {
public:
    bool intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const override;
};

} // namespace rt


