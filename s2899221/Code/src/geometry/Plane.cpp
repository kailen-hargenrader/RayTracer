#include "raytracer/geometry/Plane.h"
#include <algorithm>
#include <array>
#include <cmath>

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

    // Robust inside test for arbitrary corner ordering:
    // 1) Project to 2D based on dominant normal axis
    // 2) Order the four corners around their centroid
    // 3) Perform a convex quad edge-sign test
    auto project2D = [&](const Vec3& v)->std::pair<float,float>{
        Vec3 nabs = { std::abs(normal.x), std::abs(normal.y), std::abs(normal.z) };
        if (nabs.z >= nabs.x && nabs.z >= nabs.y) { // project to XY
            return { v.x, v.y };
        } else if (nabs.y >= nabs.x) { // XZ
            return { v.x, v.z };
        } else { // YZ
            return { v.y, v.z };
        }
    };

    std::array<std::pair<float,float>,4> c2d;
    for (int i = 0; i < 4; ++i) c2d[i] = project2D(m_corners[i]);
    auto p2d = project2D(pHit);

    // Centroid in 2D
    float cx = 0.0f, cy = 0.0f;
    for (int i = 0; i < 4; ++i) { cx += c2d[i].first; cy += c2d[i].second; }
    cx *= 0.25f; cy *= 0.25f;

    // Indices ordered by angle around centroid
    std::array<int,4> idx {0,1,2,3};
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
        float angA = std::atan2(c2d[a].second - cy, c2d[a].first - cx);
        float angB = std::atan2(c2d[b].second - cy, c2d[b].first - cx);
        return angA < angB;
    });

    std::array<Vec3,4> orderedCorners {
        m_corners[idx[0]], m_corners[idx[1]], m_corners[idx[2]], m_corners[idx[3]]
    };

    bool inside = pointInQuad(pHit, orderedCorners, normal);
    if (!inside) return false;

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


