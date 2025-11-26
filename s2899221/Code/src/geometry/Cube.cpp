#include "raytracer/geometry/Cube.h"
#include <cmath>

namespace rt {

static inline float sgn(float v) { return v < 0.0f ? -1.0f : 1.0f; }

bool Cube::intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const {
    // Object space ray
    Ray r;
    r.origin = m_worldToObject.transformPoint(rayWorld.origin);
    r.direction = m_worldToObject.transformDirection(rayWorld.direction);

    // Slab method for bounds [-0.5, 0.5]
    const float minB[3] = {-0.5f, -0.5f, -0.5f};
    const float maxB[3] = { 0.5f,  0.5f,  0.5f};
    float t0 = 0.0f;
    float t1 = 1e30f;
    int hitAxis = -1;
    int hitSign = 1;

    for (int axis = 0; axis < 3; ++axis) {
        float ro = r.origin[axis];
        float rd = r.direction[axis];
        float inv = (std::abs(rd) > 1e-8f) ? 1.0f / rd : 1e30f;
        float tNear = (minB[axis] - ro) * inv;
        float tFar  = (maxB[axis] - ro) * inv;
        int sign = 1;
        if (tNear > tFar) { std::swap(tNear, tFar); sign = -1; }
        if (tNear > t0) { t0 = tNear; hitAxis = axis; hitSign = sign; }
        if (tFar < t1) { t1 = tFar; }
        if (t0 > t1) return false;
    }

    auto computeHit = [&](float tObj, int axis)->bool{
        if (tObj <= 0.0f) return false;
        Vec3 pObj = r.at(tObj);
        Vec3 pWorld = m_objectToWorld.transformPoint(pObj);
        float tWorld = Vec3::dot(pWorld - rayWorld.origin, rayWorld.direction) / Vec3::dot(rayWorld.direction, rayWorld.direction);
        if (tWorld < tMin || tWorld > tMax) return false;
        Vec3 nObj(0.0f, 0.0f, 0.0f);
        nObj[axis] = static_cast<float>(hitSign);
        Vec3 nWorld = m_normalToWorld.transformNormal(nObj).normalized();
        hit.didHit = true;
        hit.t = tWorld;
        hit.ray = rayWorld;
        hit.position = pWorld;
        hit.setFaceNormal(nWorld);
        // Basic planar UVs per face (object space)
        if (axis == 0) { // X face
            hit.uv = { (pObj.z + 0.5f), (pObj.y + 0.5f) };
        } else if (axis == 1) { // Y face
            hit.uv = { (pObj.x + 0.5f), (pObj.z + 0.5f) };
        } else { // Z face
            hit.uv = { (pObj.x + 0.5f), (pObj.y + 0.5f) };
        }
        return true;
    };

    // Prefer nearest
    if (computeHit(t0, hitAxis)) return true;
    // Fallback to exit point
    return computeHit(t1, hitAxis);
}

} // namespace rt


