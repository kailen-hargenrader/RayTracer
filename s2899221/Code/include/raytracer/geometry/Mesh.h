#pragma once

#include <memory>
#include "raytracer/math/Mat4.h"
#include "raytracer/utils/Ray.h"
#include "raytracer/utils/Hit.h"

namespace rt {

class Mesh {
public:
    virtual ~Mesh() = default;

    // Intersect ray with primitive; if hit within [tMin, tMax], update 'hit' and return true.
    virtual bool intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const = 0;

    // Transform handling
    void setTransform(const Mat4& objectToWorld);
    void setTRS(const Vec3& translation, const Vec3& eulerXYZRadians, const Vec3& scale);

protected:
    Mat4 m_objectToWorld { Mat4::identity() };
    Mat4 m_worldToObject { Mat4::identity() };
    Mat4 m_normalToWorld { Mat4::identity() }; // (inverse transpose) for normals

    // Helpers
    void updateDerivedTransforms();
};

using MeshPtr = std::shared_ptr<Mesh>;

} // namespace rt

// Inline implementations
inline void rt::Mesh::setTransform(const Mat4& objectToWorld) {
    m_objectToWorld = objectToWorld;
    updateDerivedTransforms();
}

inline void rt::Mesh::setTRS(const Vec3& translation, const Vec3& eulerXYZRadians, const Vec3& scale) {
    m_objectToWorld = Mat4::TRS(translation, eulerXYZRadians, scale);
    updateDerivedTransforms();
}

inline void rt::Mesh::updateDerivedTransforms() {
    m_worldToObject = m_objectToWorld.inverseAffine();
    // For normals we want (inverse(M)).T acting on normals; Mat4::transformNormal will treat
    // this matrix as if it is inverse(M) and multiply by its transpose effectively.
    m_normalToWorld = m_worldToObject; // used through transformNormal
}


