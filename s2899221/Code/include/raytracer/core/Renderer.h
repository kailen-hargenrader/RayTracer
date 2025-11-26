#pragma once

#include <memory>
#include <string>
#include "raytracer/core/Scene.h"
#include "raytracer/utils/Image.h"

namespace rt {

struct RendererOptions {
    bool useBVH {true};
    int samplesPerPixel {1};
    int maxDepth {4};
    float minThroughput {0.01f};
    Vec3 ambientColor {0.08f, 0.08f, 0.08f};
    float lightIntensityScale {50.0f};
};

class Renderer {
public:
    explicit Renderer(const RendererOptions& opts) : m_opts(opts) {}

    void setScene(std::shared_ptr<Scene> scene);
    bool renderToPPM(const std::string& outputPath);

private:
    std::shared_ptr<Scene> m_scene;
    RendererOptions m_opts;

    bool intersectScene(const Ray& ray, float tMin, float tMax, Hit& hit, SceneObject const** outObj) const;
    bool isOccluded(const Vec3& origin, const Vec3& toLight, float maxDist) const;

    Vec3 shadeLocal(const Hit& hit, const SceneObject& obj, const Vec3& viewDir) const;
    Vec3 traceRay(const Ray& ray, int depth, float throughput) const;

    static Vec3 reflect(const Vec3& v, const Vec3& n) {
        return v - 2.0f * Vec3::dot(v, n) * n;
    }
    static bool refract(const Vec3& v, const Vec3& n, float eta, Vec3& out) {
        float cosi = -std::max(-1.0f, std::min(1.0f, Vec3::dot(v, n)));
        float sint2 = eta*eta * (1.0f - cosi*cosi);
        if (sint2 > 1.0f) return false;
        float cost = std::sqrt(std::max(0.0f, 1.0f - sint2));
        out = eta * v + (eta * cosi - cost) * n;
        return true;
    }
    static float fresnelSchlick(float cosTheta, float F0) {
        return F0 + (1.0f - F0) * std::pow(1.0f - cosTheta, 5.0f);
    }
    static float shininessFromRoughness(float r) {
        // Map [0,1] to [2,256]
        r = std::max(0.0f, std::min(1.0f, r));
        return 2.0f + (1.0f - r) * 254.0f;
    }
};

} // namespace rt


