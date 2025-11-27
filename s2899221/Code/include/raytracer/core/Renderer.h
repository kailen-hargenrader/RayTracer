#pragma once

#include <memory>
#include <vector>
#include <random>
#include "raytracer/core/Scene.h"
#include "raytracer/utils/Image.h"
#include "raytracer/utils/Ray.h"

namespace rt {

struct RenderOptions {
	bool useBVH { true };
	int samplesPerPixel { 1 };
	int maxDepth { 5 }; // maximum path length in segments; 1 = primary only (no bounces)
	// When true, disables both reflection and refraction. Transparent materials pass rays straight through
	// (no bending), and opaque materials render with local lighting only (no glossy reflection).
	bool transmissionOnly { false };
	// Maximum number of stochastic reflection samples used to blur reflections based on material roughness.
	// 1 = no blur (perfect mirror). Higher values increase quality and cost.
	int roughnessMaxSamples { 1 };
	bool useRussianRoulette { false }; // reserved
	float minThroughput { 0.02f }; // reserved
	Vec3 backgroundColor { 0.04f, 0.04f, 0.02f };
};

class Renderer {
public:
	void render(const Scene& scene, Image& outImage, const RenderOptions& opts);

private:
	struct TraceResult {
		bool hit { false };
		Vec3 color { 0.0f, 0.0f, 0.0f };
	};

	TraceResult traceRay(const Scene& scene, const Ray& ray, const RenderOptions& opts, int depth, std::mt19937& rng) const;
	Vec3 shadeBlinnPhong(const Scene& scene, const Hit& hit, const Material& material, const RenderOptions& opts) const;
	Vec3 backgroundColor(const Ray& ray, const RenderOptions& opts) const;
	// Returns [0..1] multiplicative light visibility along dirToLight, accounting for transparency
	float shadowTransmittance(const Scene& scene, const Vec3& pos, const Vec3& dirToLight, float maxDist) const;
};

} // namespace rt


