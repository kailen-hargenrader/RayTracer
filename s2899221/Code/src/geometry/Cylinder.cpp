#include "raytracer/geometry/Cylinder.h"
#include <cmath>

namespace rt {

bool Cylinder::intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const {
    // Transform ray to object space
    Ray r;
    r.origin = m_worldToObject.transformPoint(rayWorld.origin);
    r.direction = m_worldToObject.transformDirection(rayWorld.direction);

    const float radius = 0.5f;
    const float yMin = -0.5f, yMax = 0.5f;

    bool anyHit = false;
    float bestTWorld = tMax;
    Hit bestHit = hit;

    auto commit = [&](float tWorld, const Vec3& pWorld, const Vec3& outwardNormalWorld, const Vec2& uv)->void{
        if (tWorld < bestTWorld && tWorld >= tMin && tWorld <= tMax) {
            anyHit = true;
            bestTWorld = tWorld;
            bestHit.didHit = true;
            bestHit.t = tWorld;
            bestHit.ray = rayWorld;
            bestHit.position = pWorld;
            bestHit.setFaceNormal(outwardNormalWorld);
            bestHit.uv = uv;
        }
    };

    // Side hit: x^2 + z^2 = r^2, y within limits
    {
        const float dx = r.direction.x;
        const float dz = r.direction.z;
        const float ox = r.origin.x;
        const float oz = r.origin.z;
        float a = dx*dx + dz*dz;
        float b = 2.0f * (ox*dx + oz*dz);
        float c = ox*ox + oz*oz - radius*radius;
        float disc = b*b - 4.0f*a*c;
        if (disc >= 0.0f && a > 1e-8f) {
            float sqrtD = std::sqrt(disc);
            float t0 = (-b - sqrtD) / (2.0f * a);
            float t1 = (-b + sqrtD) / (2.0f * a);
            auto testSideRoot = [&](float tObj){
                if (tObj <= 0.0f) return;
                Vec3 pObj = r.at(tObj);
                if (pObj.y < yMin - 1e-6f || pObj.y > yMax + 1e-6f) return;
                Vec3 pWorld = m_objectToWorld.transformPoint(pObj);
                float tWorld = Vec3::dot(pWorld - rayWorld.origin, rayWorld.direction) / Vec3::dot(rayWorld.direction, rayWorld.direction);
                Vec3 nObj = Vec3(pObj.x, 0.0f, pObj.z).normalized();
                Vec3 nWorld = m_normalToWorld.transformNormal(nObj).normalized();
                // UV around side
                float u = 0.5f + std::atan2(nObj.z, nObj.x) / (2.0f * 3.1415926535f);
                float v = (pObj.y - yMin) / (yMax - yMin);
                commit(tWorld, pWorld, nWorld, {u, v});
            };
            testSideRoot(t0);
            testSideRoot(t1);
        }
    }

    // Caps: y = yMin and y = yMax, circle x^2 + z^2 <= r^2
    auto testCap = [&](float yPlane, const Vec3& nObj){
        float denom = r.direction.y;
        if (std::abs(denom) < 1e-8f) return;
        float tObj = (yPlane - r.origin.y) / denom;
        if (tObj <= 0.0f) return;
        Vec3 pObj = r.at(tObj);
        if (pObj.x*pObj.x + pObj.z*pObj.z > radius*radius + 1e-6f) return;
        Vec3 pWorld = m_objectToWorld.transformPoint(pObj);
        float tWorld = Vec3::dot(pWorld - rayWorld.origin, rayWorld.direction) / Vec3::dot(rayWorld.direction, rayWorld.direction);
        Vec3 nWorld = m_normalToWorld.transformNormal(nObj).normalized();
        // UV for cap (mapped into [0,1])
        float u = (pObj.x + radius) / (2.0f * radius);
        float v = (pObj.z + radius) / (2.0f * radius);
        commit(tWorld, pWorld, nWorld, {u, v});
    };
    testCap(yMin, {0.0f, -1.0f, 0.0f});
    testCap(yMax, {0.0f,  1.0f, 0.0f});

    if (anyHit) {
        hit = bestHit;
        return true;
    }
    return false;
}

} // namespace rt


