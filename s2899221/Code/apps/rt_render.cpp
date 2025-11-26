#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include "raytracer/core/Scene.h"
#include "raytracer/core/Renderer.h"

using namespace rt;

static void printUsage() {
    std::cout << "Usage:\n"
              << "  rt_render <scene.json> <output.ppm> [--bvh|--no-bvh] [--spp N] [--maxDepth N] [--minThroughput X] [--ambient r,g,b] [--lightScale X]\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printUsage();
        return 1;
    }
    std::string jsonPath = argv[1];
    std::string outPath = argv[2];

    RendererOptions opts;
    // Defaults
    opts.useBVH = true;
    opts.samplesPerPixel = 4;
    opts.maxDepth = 5;
    opts.minThroughput = 0.02f;

    for (int i = 3; i < argc; ++i) {
        if (std::strcmp(argv[i], "--bvh") == 0) {
            opts.useBVH = true;
        } else if (std::strcmp(argv[i], "--no-bvh") == 0) {
            opts.useBVH = false;
        } else if (std::strcmp(argv[i], "--spp") == 0 && i + 1 < argc) {
            opts.samplesPerPixel = std::max(1, std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--maxDepth") == 0 && i + 1 < argc) {
            opts.maxDepth = std::max(0, std::atoi(argv[++i]));
        } else if (std::strcmp(argv[i], "--minThroughput") == 0 && i + 1 < argc) {
            opts.minThroughput = std::max(0.0f, static_cast<float>(std::atof(argv[++i])));
        } else if (std::strcmp(argv[i], "--ambient") == 0 && i + 1 < argc) {
            std::string s = argv[++i];
            float r=0,g=0,b=0;
            if (std::sscanf(s.c_str(), "%f,%f,%f", &r, &g, &b) == 3) {
                opts.ambientColor = {r,g,b};
            }
        } else if (std::strcmp(argv[i], "--lightScale") == 0 && i + 1 < argc) {
            opts.lightIntensityScale = static_cast<float>(std::atof(argv[++i]));
        } else {
            std::cout << "Unknown or incomplete option: " << argv[i] << "\n";
            printUsage();
            return 1;
        }
    }

    auto scene = Scene::loadFromJsonFile(jsonPath);
    if (!scene) {
        std::cerr << "Failed to load scene: " << jsonPath << "\n";
        return 2;
    }

    Renderer renderer(opts);
    renderer.setScene(scene);

    std::filesystem::create_directories(std::filesystem::path(outPath).parent_path());

    if (!renderer.renderToPPM(outPath)) {
        std::cerr << "Render failed.\n";
        return 3;
    }
    std::cout << "Wrote: " << outPath << "\n";
    return 0;
}


