#include "raytracer/geometry/Sphere.h"

namespace rt {

bool Sphere::intersect(const Ray& rayWorld, float tMin, float tMax, Hit& hit) const {
    // Transform ray to object space
    Ray rayObj;
    rayObj.origin = m_worldToObject.transformPoint(rayWorld.origin);
    rayObj.direction = m_worldToObject.transformDirection(rayWorld.direction);

    // Sphere centered at origin radius 0.5
    const float r = 0.5f;
    const float a = Vec3::dot(rayObj.direction, rayObj.direction);
    const Vec3 oc = rayObj.origin; // since center is (0,0,0)
    const float b = 2.0f * Vec3::dot(oc, rayObj.direction);
    const float c = Vec3::dot(oc, oc) - r*r;
    const float disc = b*b - 4.0f*a*c;
    if (disc < 0.0f) return false;
    const float sqrtD = std::sqrt(disc);

    auto tryRoot = [&](float tObj)->bool{
        if (tObj <= 0.0f) return false;
        // Convert to world-space t via projection onto world ray direction
        Vec3 pObj = rayObj.at(tObj);
        Vec3 pWorld = m_objectToWorld.transformPoint(pObj);
        float tWorld = Vec3::dot(pWorld - rayWorld.origin, rayWorld.direction) / Vec3::dot(rayWorld.direction, rayWorld.direction);
        if (tWorld < tMin || tWorld > tMax) return false;

        // Fill hit
        Vec3 outwardNormalObj = (pObj * (1.0f / r)).normalized();
        Vec3 outwardNormalWorld = m_normalToWorld.transformNormal(outwardNormalObj).normalized();

        hit.didHit = true;
        hit.t = tWorld;
        hit.ray = rayWorld;
        hit.position = pWorld;
        hit.setFaceNormal(outwardNormalWorld);
        // UV for sphere: Blender-style (v increases upward)
        float u = 0.5f + std::atan2(outwardNormalObj.z, outwardNormalObj.x) / (2.0f * 3.1415926535f);
        float v = 0.5f + std::asin(outwardNormalObj.y) / 3.1415926535f;
        hit.uv = {u, v};
        return true;
    };

    float root = (-b - sqrtD) / (2.0f * a);
    if (tryRoot(root)) return true;
    root = (-b + sqrtD) / (2.0f * a);
    return tryRoot(root);
}

} // namespace rt


