#pragma once

#include "raytracer/geometry/Mesh.h"

namespace rt {

// Canonical axis-aligned cube centered at origin with side length 1 (bounds [-0.5,0.5])
class Cube : public Mesh {
public:
    bool intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const override;
};

} // namespace rt


