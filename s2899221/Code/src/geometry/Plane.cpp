#include "raytracer/geometry/Plane.h"

namespace rt {

static inline bool pointInQuad(const Vec3& p, const std::array<Vec3,4>& c, const Vec3& n) {
    // All edge cross-products should have same sign relative to normal
    for (int i = 0; i < 4; ++i) {
        const Vec3& a = c[i];
        const Vec3& b = c[(i+1)&3];
        Vec3 edge = b - a;
        Vec3 toP = p - a;
        Vec3 cprod = Vec3::cross(edge, toP);
        if (Vec3::dot(cprod, n) < -1e-6f) return false;
    }
    return true;
}

bool Plane::intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const {
    if (m_corners[0].lengthSquared() == 0.0f && m_corners[1].lengthSquared() == 0.0f &&
        m_corners[2].lengthSquared() == 0.0f && m_corners[3].lengthSquared() == 0.0f) {
        return false; // not initialized
    }
    const Vec3& p0 = m_corners[0];
    Vec3 e1 = m_corners[1] - p0;
    Vec3 e2 = m_corners[2] - p0;
    Vec3 normal = Vec3::cross(e1, e2).normalized();
    float denom = Vec3::dot(rayWorld.direction, normal);
    if (std::abs(denom) < 1e-8f) return false; // parallel
    float t = Vec3::dot(p0 - rayWorld.origin, normal) / denom;
    if (t < tMin || t > tMax) return false;
    Vec3 pHit = rayWorld.at(t);
    if (!pointInQuad(pHit, m_corners, normal)) return false;

    hit.didHit = true;
    hit.t = t;
    hit.ray = rayWorld;
    hit.position = pHit;
    hit.setFaceNormal(normal);

    // Simple planar UV via projection (choose dominant axis of normal)
    Vec3 n = normal;
    n = {std::abs(n.x), std::abs(n.y), std::abs(n.z)};
    Vec3 pLocal = pHit - p0;
    float u=0.0f, v=0.0f;
    if (n.z >= n.x && n.z >= n.y) { // project to XY
        u = (pLocal.x - 0.0f);
        v = (pLocal.y - 0.0f);
    } else if (n.y >= n.x) { // XZ
        u = (pLocal.x - 0.0f);
        v = (pLocal.z - 0.0f);
    } else { // YZ
        u = (pLocal.y - 0.0f);
        v = (pLocal.z - 0.0f);
    }
    hit.uv = {u, v};
    return true;
}

} // namespace rt


