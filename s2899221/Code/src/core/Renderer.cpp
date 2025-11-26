#include "raytracer/core/Renderer.h"
#include "raytracer/utils/Texture.h"
#include <random>

namespace rt {

void Renderer::setScene(std::shared_ptr<Scene> scene) {
    m_scene = std::move(scene);
    if (m_scene && m_opts.useBVH) {
        m_scene->bvh.build(m_scene->objects);
    }
    // Load textures if any
    if (m_scene) {
        for (auto& o : m_scene->objects) {
            if (!o.material.baseColorTexturePath.empty()) {
                o.material.baseColorTexture = std::make_shared<Texture>();
                o.material.baseColorTexture->loadFromFile(o.material.baseColorTexturePath);
            }
        }
    }
}

bool Renderer::intersectScene(const Ray& ray, float tMin, float tMax, Hit& hit, SceneObject const** outObj) const {
    bool any = false;
    if (m_opts.useBVH) {
        int objIndex = -1;
        any = m_scene->bvh.intersect(m_scene->objects, ray, tMin, tMax, hit, &objIndex);
        if (any && outObj) {
            *outObj = (objIndex >= 0 && objIndex < static_cast<int>(m_scene->objects.size()))
                        ? &m_scene->objects[objIndex]
                        : nullptr;
        }
        return any;
    }
    // Unaccelerated
    float closest = tMax;
    const SceneObject* found = nullptr;
    Hit tmp;
    for (const auto& obj : m_scene->objects) {
        if (!obj.bounds.intersect(ray, tMin, closest)) continue;
        if (obj.mesh->intersect(ray, tMin, closest, tmp)) {
            any = true;
            closest = tmp.t;
            hit = tmp;
            found = &obj;
        }
    }
    if (outObj) *outObj = found;
    return any;
}

bool Renderer::isOccluded(const Vec3& origin, const Vec3& toLight, float maxDist) const {
    Vec3 dir = toLight.normalized();
    Ray shadowRay(origin + dir * 1e-4f, dir);
    Hit h;
    SceneObject const* obj = nullptr;
    if (m_opts.useBVH) {
        return m_scene->bvh.intersect(m_scene->objects, shadowRay, 1e-4f, maxDist - 1e-4f, h);
    } else {
        for (const auto& o : m_scene->objects) {
            if (!o.bounds.intersect(shadowRay, 1e-4f, maxDist - 1e-4f)) continue;
            Hit tmp;
            if (o.mesh->intersect(shadowRay, 1e-4f, maxDist - 1e-4f, tmp)) {
                return true;
            }
        }
        return false;
    }
}

Vec3 Renderer::shadeLocal(const Hit& hit, const SceneObject& obj, const Vec3& viewDir) const {
    const Material& m = obj.material;
    // Sample baseColor texture if present
    Vec3 base = m.baseColor;
    if (m.baseColorTexture && hit.uv.x >= 0.0f) {
        base = m.baseColorTexture->sample(hit.uv.x, hit.uv.y);
    }
    // Disney-ish mixing
    float metallic = std::max(0.0f, std::min(1.0f, m.metallic));
    Vec3 F0Color = (1.0f - metallic) * Vec3(0.04f, 0.04f, 0.04f) + metallic * base;
    Vec3 kd = base * (1.0f - metallic);
    float shininess = shininessFromRoughness(m.roughness);
    Vec3 color = m_opts.ambientColor * base;
    for (const auto& L : m_scene->lights) {
        Vec3 toL = L.position - hit.position;
        float dist2 = std::max(1e-4f, Vec3::dot(toL, toL));
        float dist = std::sqrt(dist2);
        if (isOccluded(hit.position, toL, dist)) continue;
        Vec3 ldir = toL / dist;
        float NdotL = std::max(0.0f, Vec3::dot(hit.normal, ldir));
        Vec3 diffuse = kd * NdotL;
        Vec3 hdir = (ldir + viewDir).normalized();
        float NdotH = std::max(0.0f, Vec3::dot(hit.normal, hdir));
        Vec3 specular = F0Color * std::pow(NdotH, shininess);
        float attenuation = (L.intensity * m_opts.lightIntensityScale) / (dist2);
        color += (diffuse + specular) * attenuation;
    }
    return Vec3::clamp01(color);
}

Vec3 Renderer::traceRay(const Ray& ray, int depth, float throughput) const {
    if (depth > m_opts.maxDepth || throughput < m_opts.minThroughput) return {0,0,0};

    Hit hit;
    const SceneObject* obj = nullptr;
    if (!intersectScene(ray, 1e-4f, 1e30f, hit, &obj) || obj == nullptr) {
        // Simple sky gradient background (linear)
        float t = 0.5f * (ray.direction.z + 1.0f);
        Vec3 cTop = {0.6f, 0.8f, 1.0f};
        Vec3 cBot = {1.0f, 1.0f, 1.0f};
        return (1.0f - t) * cBot + t * cTop;
    }
    Vec3 viewDir = (-ray.direction).normalized();

    // Local shading
    Vec3 localColor = shadeLocal(hit, *obj, viewDir);

    // Fresnel and reflection/refraction
    const Material& mat = obj->material;
    float ior = std::max(1.0f, mat.ior);
    // F0 for dielectrics from ior; for metals, approximate using base color luminance
    float F0 = (mat.metallic > 0.5f)
        ? std::max(0.02f, (mat.baseColor.x + mat.baseColor.y + mat.baseColor.z) / 3.0f)
        : ((1.0f - ior) / (1.0f + ior)); // We'll square later
    F0 = F0 * F0;

    float cosTheta = std::max(0.0f, Vec3::dot(hit.normal, viewDir));
    float Fr = fresnelSchlick(cosTheta, F0);

    Vec3 reflectedColor(0,0,0), refractedColor(0,0,0);
    // Reflection
    Vec3 reflDir = reflect(-viewDir, hit.normal).normalized();
    reflectedColor = traceRay(Ray(hit.position + reflDir * 1e-4f, reflDir), depth + 1, throughput * Fr);

    // Refraction only if not fully opaque
    if (mat.alpha < 1.0f) {
        Vec3 n = hit.normal;
        float etai = 1.0f, etat = ior;
        float cosi = std::max(-1.0f, std::min(1.0f, Vec3::dot(viewDir, n)));
        if (cosi < 0.0f) {
            cosi = -cosi;
        } else {
            std::swap(etai, etat);
            n = -n;
        }
        float eta = etai / etat;
        Vec3 refrDir;
        if (refract(-viewDir, n, eta, refrDir)) {
            refrDir = refrDir.normalized();
            refractedColor = traceRay(Ray(hit.position + refrDir * 1e-4f, refrDir), depth + 1, throughput * (1.0f - Fr));
        }
    }

    // Energy-conserving mix: dielectric uses Fr; metallic reduces diffuse already
    Vec3 result = localColor * (1.0f - Fr) + reflectedColor * Fr + refractedColor * (1.0f - Fr) * (1.0f - mat.alpha);
    return Vec3::clamp01(result);
}

bool Renderer::renderToPPM(const std::string& outputPath) {
    if (!m_scene) return false;
    int W = m_scene->camera.width();
    int H = m_scene->camera.height();
    Image img(W, H);

    std::mt19937 rng(12345);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    int spp = std::max(1, m_opts.samplesPerPixel);

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Vec3 color(0,0,0);
            for (int s = 0; s < spp; ++s) {
                float jx = (spp == 1) ? 0.5f : dist(rng);
                float jy = (spp == 1) ? 0.5f : dist(rng);
                Ray ray = m_scene->camera.generateRay(x, y, jx, jy);
                color += traceRay(ray, 0, 1.0f);
            }
            color /= static_cast<float>(spp);
            // Gamma encode to approximate sRGB
            color = { std::pow(color.x, 1.0f/2.2f), std::pow(color.y, 1.0f/2.2f), std::pow(color.z, 1.0f/2.2f) };
            img.setPixel(x, y, color);
        }
    }
    return img.writePPM(outputPath, true);
}

} // namespace rt


