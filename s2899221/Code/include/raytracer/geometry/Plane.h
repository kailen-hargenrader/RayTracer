#pragma once

#include "raytracer/geometry/Mesh.h"
#include <array>

namespace rt {

// Plane defined by four world-space corners (quad), in order
class Plane : public Mesh {
public:
    Plane() = default;
    explicit Plane(const std::array<Vec3,4>& cornersWorld) : m_corners(cornersWorld) {}

    void setCorners(const std::array<Vec3,4>& cornersWorld) { m_corners = cornersWorld; }
    const std::array<Vec3,4>& corners() const { return m_corners; }

    bool intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const override;

private:
    std::array<Vec3,4> m_corners; // world-space quad
};

} // namespace rt


