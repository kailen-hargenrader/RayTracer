#pragma once

#include "raytracer/geometry/Mesh.h"

namespace rt {

// Canonical finite cylinder along Y axis: radius 0.5, y in [-0.5, 0.5]
class Cylinder : public Mesh {
public:
    bool intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const override;
};

} // namespace rt


