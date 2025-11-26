#include <iostream>
#include <string>
#include <memory>
#include <cstdlib>
#include <algorithm>
#include "raytracer/core/Scene.h"
#include "raytracer/core/Renderer.h"
#include "raytracer/utils/Image.h"

using namespace rt;

static void printUsage() {
	std::cout << "Usage: rt_render <scene.json> <out.ppm> [--bvh] [--spp N] [--maxDepth N] [--minThroughput f] [--no-rr]\n";
}

int main(int argc, char** argv) {
	if (argc < 3) {
		printUsage();
		return 1;
	}
	std::string scenePath = argv[1];
	std::string outPath = argv[2];

	RenderOptions opts;
	// Defaults aligned with batch
	opts.useBVH = false;
	opts.samplesPerPixel = 1;
	opts.maxDepth = 1;
	opts.useRussianRoulette = false;
	opts.minThroughput = 0.02f;

	for (int i = 3; i < argc; ++i) {
		std::string a = argv[i];
		if (a == "--bvh") {
			opts.useBVH = true;
		} else if (a == "--no-rr") {
			opts.useRussianRoulette = false;
		} else if (a == "--spp" && i + 1 < argc) {
			opts.samplesPerPixel = std::max(1, std::atoi(argv[++i]));
		} else if (a == "--maxDepth" && i + 1 < argc) {
			opts.maxDepth = std::max(1, std::atoi(argv[++i]));
		} else if (a == "--minThroughput" && i + 1 < argc) {
			opts.minThroughput = std::max(0.0f, static_cast<float>(std::atof(argv[++i])));
		} else {
			std::cout << "Unknown or malformed option: " << a << "\n";
			printUsage();
			return 1;
		}
	}

	auto scene = Scene::loadFromJsonFile(scenePath);
	if (!scene) {
		std::cerr << "Failed to load scene: " << scenePath << "\n";
		return 1;
	}

	Image img(scene->camera.width(), scene->camera.height());
	Renderer renderer;
	renderer.render(*scene, img, opts);
	if (!img.writePPM(outPath, true)) {
		std::cerr << "Failed to write image: " << outPath << "\n";
		return 1;
	}

	std::cout << "Rendered " << scene->camera.width() << "x" << scene->camera.height()
		<< " spp=" << opts.samplesPerPixel << " to " << outPath << "\n";
	return 0;
}


