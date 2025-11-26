#include "raytracer/core/Renderer.h"
#include "raytracer/utils/Texture.h"
#include <random>
#include <limits>
#include <iostream>
#include <algorithm>

namespace rt {

static constexpr float kEpsilon = 1e-4f;

static inline Vec3 reflect(const Vec3& v, const Vec3& n) {
	return v - n * (2.0f * Vec3::dot(v, n));
}

static inline bool refract(const Vec3& v, const Vec3& n, float eta, Vec3& refractedOut) {
	// v and n should be normalized; eta = eta_i/eta_t
	float cosTheta = std::max(0.0f, -Vec3::dot(v, n));
	Vec3 rOutPerp = (v + n * cosTheta) * eta;
	float k = 1.0f - rOutPerp.lengthSquared();
	if (k < 0.0f) {
		return false; // total internal reflection
	}
	Vec3 rOutParallel = n * -std::sqrt(k);
	refractedOut = rOutPerp + rOutParallel;
	return true;
}

static inline float schlickReflectance(float cosTheta, float ior) {
	// cosTheta is cosine of angle between incident and normal (assumed in [0,1])
	float r0 = (1.0f - ior) / (1.0f + ior);
	r0 = r0 * r0;
	return r0 + (1.0f - r0) * std::pow(1.0f - cosTheta, 5.0f);
}

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

	// Reflection / Refraction
	Vec3 finalColor = local;
	if (depth < opts.maxDepth) {
		const Vec3 N = hit.normal.normalized();
		const Vec3 I = hit.ray.direction.normalized();

		const float metallic = saturate(mat.metallic);
		const float alpha = saturate(mat.alpha);
		const float transparency = 1.0f - alpha;
		const float ior = std::max(1.0f, mat.ior);

		// Fresnel term (dielectric); for metals we bias toward reflection
		const float cosTheta = std::max(0.0f, -Vec3::dot(I, N));
		float F = schlickReflectance(cosTheta, ior);
		F = F * (1.0f - metallic) + 1.0f * metallic; // metals reflect nearly all

		// Opaque: mix local with reflection
		if (transparency <= 1e-3f) {
			Vec3 Rdir = reflect(I, N).normalized();
			Ray rRay(hit.position + N * kEpsilon * 4.0f, Rdir);
			TraceResult rTr = traceRay(scene, rRay, opts, depth + 1);
			finalColor = (1.0f - F) * local + F * rTr.color;
		} else {
			// Simplest transparency: continue the ray straight through without bending
			Vec3 passDir = I;
			Ray passRay(hit.position + passDir * kEpsilon * 8.0f, passDir);
			Vec3 passCol = traceRay(scene, passRay, opts, depth + 1).color;
			finalColor = alpha * local + transparency * passCol;
		}
	}

	result.hit = true;
	result.color = saturate(finalColor);
	return result;
}

Vec3 Renderer::shadeBlinnPhong(const Scene& scene, const Hit& hit, const Material& material, const RenderOptions& opts) const {
	// Basic parameters
	const Vec3 N = hit.normal.normalized();
	const Vec3 V = (-hit.ray.direction).normalized();

	// Base color: sample texture if available, otherwise use constant
	Vec3 baseColor = material.baseColor;
	if (material.baseColorTexture && material.baseColorTexture->isValid()) {
		// Blender-style UVs appear swapped relative to our sampling; use (v, u)
		baseColor = material.baseColorTexture->sample(hit.uv.y, hit.uv.x);
	}
	// Diffuse color and specular color based on metallic
	const float metallic = saturate(material.metallic);
	const float roughness = saturate(material.roughness);

	// Map roughness to shininess exponent (simple heuristic)
	// roughness=0 -> high exponent, roughness=1 -> low exponent
	const float gloss = (1.0f - roughness);
	const float shininess = 2.0f + gloss * gloss * 256.0f;

	// Simple specular color model (keep metallic workflow, but normalize Blinn-Phong)
	const Vec3 F0_dielectric(0.04f, 0.04f, 0.04f);
	const Vec3 specularColor = F0_dielectric * (1.0f - metallic) + baseColor * metallic;
	const Vec3 diffuseColor = baseColor * (1.0f - metallic);

	// Normalized Lambert + Blinn-Phong, slightly higher ambient for Blender-like fill
	const float invPi = 0.31830988618f; // 1/pi
	const float specNorm = (shininess + 8.0f) * (1.0f / (8.0f * 3.1415926535f));
	const float ambientK = 0.08f;

	Vec3 color = diffuseColor * ambientK;

	for (const PointLight& Ls : scene.lights) {
		Vec3 L = (Ls.position - hit.position);
		const float dist2 = std::max(1e-6f, L.lengthSquared());
		const float dist = std::sqrt(dist2);
		L = L / dist;

		// Visibility accounting for transparency
		float vis = shadowTransmittance(scene, hit.position + N * kEpsilon * 4.0f, L, dist - kEpsilon * 8.0f);
		if (vis <= 1e-4f) continue;

		const float NdotL = std::max(0.0f, Vec3::dot(N, L));
		Vec3 H = (V + L).normalized();
		const float NdotH = std::max(0.0f, Vec3::dot(N, H));

		// Normalized Blinn-Phong terms
		Vec3 diffuse = diffuseColor * (NdotL * invPi);
		Vec3 specular = specularColor * (specNorm * std::pow(NdotH, shininess) * NdotL);

		// Simple inverse-square attenuation using radiant intensity proxy
		const float atten = Ls.intensity / (4.0f * 3.1415926535f * dist2);
		color += (diffuse + specular) * atten * vis;
	}

	return color;
}

Vec3 Renderer::backgroundColor(const Ray& ray, const RenderOptions& opts) const {
	// Uniform background color
	return opts.backgroundColor;
}

float Renderer::shadowTransmittance(const Scene& scene, const Vec3& pos, const Vec3& dirToLight, float maxDist) const {
	Ray ray(pos, dirToLight);
	float trans = 1.0f;
	float remainingDist = maxDist;

	for (int iter = 0; iter < 64 && trans > 1e-3f && remainingDist > kEpsilon; ++iter) {
		Hit h;
		int objIdx = -1;
		bool hit = scene.bvh.intersect(scene.objects, ray, kEpsilon, remainingDist, h, &objIdx);
		if (!hit || objIdx < 0) break;

		const Material& m = scene.objects[objIdx].material;
		const float alpha = saturate(m.alpha);
		const float metallic = saturate(m.metallic);

		// Fully opaque or metal blocks light
		if (alpha >= 0.999f || metallic >= 0.999f) {
			return 0.0f;
		}

		// Approximate Fresnel transmission factor for dielectrics
		float ior = std::max(1.0f, m.ior);
		float cosTheta = std::fabs(Vec3::dot(-dirToLight, h.normal.normalized()));
		float F = schlickReflectance(cosTheta, ior);
		float materialTrans = (1.0f - alpha) * (1.0f - F);

		trans *= materialTrans;
		// Advance the ray origin slightly forward and reduce remaining distance
		float advance = h.t + kEpsilon * 8.0f;
		ray.origin = ray.at(advance);
		remainingDist -= advance;
	}

	return saturate(trans);
}

} // namespace rt


