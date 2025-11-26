#include "raytracer/core/Renderer.h"
#include <random>
#include <limits>
#include <iostream>
#include <algorithm>

namespace rt {

static constexpr float kEpsilon = 1e-4f;

static inline float saturate(float x) {
	return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x);
}

static inline Vec3 saturate(const Vec3& v) {
	return { saturate(v.x), saturate(v.y), saturate(v.z) };
}

void Renderer::render(const Scene& scene, Image& outImage, const RenderOptions& opts) {
	const int W = scene.camera.width();
	const int H = scene.camera.height();
	outImage.resize(W, H);

	// Build BVH if requested
	if (opts.useBVH) {
		const_cast<Scene&>(scene).bvh.build(scene.objects);
	}

	std::mt19937 rng(42);
	std::uniform_real_distribution<float> uni(0.0f, 1.0f);

	for (int y = 0; y < H; ++y) {
		for (int x = 0; x < W; ++x) {
			Vec3 accum(0.0f, 0.0f, 0.0f);
			for (int s = 0; s < std::max(1, opts.samplesPerPixel); ++s) {
				float jx = (opts.samplesPerPixel > 1) ? uni(rng) : 0.5f;
				float jy = (opts.samplesPerPixel > 1) ? uni(rng) : 0.5f;
				Ray ray = scene.camera.generateRay(x, y, jx, jy);
				TraceResult tr = traceRay(scene, ray, opts, 0);
				accum += tr.color;
			}
			float invSpp = 1.0f / static_cast<float>(std::max(1, opts.samplesPerPixel));
			outImage.setPixel(x, y, accum * invSpp);
		}
	}
}

Renderer::TraceResult Renderer::traceRay(const Scene& scene, const Ray& ray, const RenderOptions& opts, int depth) const {
	TraceResult result;
	result.hit = false;
	result.color = backgroundColor(ray, opts);

	Hit hit;
	hit.ray = ray;
	int hitObjectIndex = -1;

	bool anyHit = false;
	if (opts.useBVH) {
		anyHit = scene.bvh.intersect(scene.objects, ray, kEpsilon, std::numeric_limits<float>::infinity(), hit, &hitObjectIndex);
	} else {
		float closest = std::numeric_limits<float>::infinity();
		for (size_t i = 0; i < scene.objects.size(); ++i) {
			Hit tmp = hit;
			if (scene.objects[i].mesh->intersect(ray, kEpsilon, closest, tmp)) {
				closest = tmp.t;
				tmp.albedo = scene.objects[i].material.baseColor;
				hit = tmp;
				hitObjectIndex = static_cast<int>(i);
				anyHit = true;
			}
		}
	}

	if (!anyHit || hitObjectIndex < 0) {
		return result;
	}

	const SceneObject& obj = scene.objects[hitObjectIndex];
	const Material& mat = obj.material;
	Vec3 local = shadeBlinnPhong(scene, hit, mat, opts);

	// Placeholders for future reflection/refraction mixing based on metallic/alpha/ior
	// For now, just output local shading
	result.hit = true;
	result.color = saturate(local);
	return result;
}

Vec3 Renderer::shadeBlinnPhong(const Scene& scene, const Hit& hit, const Material& material, const RenderOptions& opts) const {
	// Basic parameters
	const Vec3 N = hit.normal.normalized();
	const Vec3 V = (-hit.ray.direction).normalized();

	// Diffuse color and specular color based on metallic
	const Vec3 baseColor = material.baseColor;
	const float metallic = saturate(material.metallic);
	const float roughness = saturate(material.roughness);

	// Map roughness to shininess exponent (simple heuristic)
	// roughness=0 -> high exponent, roughness=1 -> low exponent
	const float gloss = (1.0f - roughness);
	const float shininess = 2.0f + gloss * gloss * 256.0f;

	const Vec3 F0_dielectric(0.04f, 0.04f, 0.04f);
	const Vec3 specularColor = F0_dielectric * (1.0f - metallic) + baseColor * metallic;
	const Vec3 diffuseColor = baseColor * (1.0f - metallic);
	const float ambientK = 0.02f;

	Vec3 color = diffuseColor * ambientK;

	for (const PointLight& Ls : scene.lights) {
		Vec3 L = (Ls.position - hit.position);
		const float dist2 = std::max(1e-6f, L.lengthSquared());
		const float dist = std::sqrt(dist2);
		L = L / dist;

		// Hard shadows
		if (occludedToLight(scene, hit.position + N * kEpsilon * 4.0f, L, dist - kEpsilon * 8.0f)) {
			continue;
		}

		const float NdotL = std::max(0.0f, Vec3::dot(N, L));
		Vec3 H = (V + L).normalized();
		const float NdotH = std::max(0.0f, Vec3::dot(N, H));

		// Blinn-Phong terms
		Vec3 diffuse = diffuseColor * NdotL;
		Vec3 specular = specularColor * std::pow(NdotH, shininess);

		// Simple inverse-square attenuation using radiant intensity proxy
		const float atten = Ls.intensity / (4.0f * 3.1415926535f * dist2);
		color += (diffuse + specular) * atten;
	}

	return color;
}

Vec3 Renderer::backgroundColor(const Ray& ray, const RenderOptions& opts) const {
	// Simple vertical gradient
	Vec3 d = ray.direction.normalized();
	float t = 0.5f * (d.z + 1.0f);
	return (1.0f - t) * opts.backgroundBottom + t * opts.backgroundTop;
}

bool Renderer::occludedToLight(const Scene& scene, const Vec3& pos, const Vec3& dirToLight, float maxDist) const {
	Ray shadowRay(pos, dirToLight);
	Hit h;
	// Traverse via BVH if available, otherwise brute-force
	if (!scene.objects.empty()) {
		// Prefer BVH (always present in Scene)
		if (scene.bvh.intersect(scene.objects, shadowRay, kEpsilon, maxDist, h, nullptr)) return true;
		// In case BVH is empty (scene not built), brute-force
		for (const SceneObject& o : scene.objects) {
			Hit tmp;
			if (o.mesh->intersect(shadowRay, kEpsilon, maxDist, tmp)) return true;
		}
	}
	return false;
}

} // namespace rt


