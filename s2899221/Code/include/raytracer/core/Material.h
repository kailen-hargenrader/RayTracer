#pragma once

#include "raytracer/math/Vec3.h"
#include <string>
#include <memory>

namespace rt {

class Texture;

struct Material {
    Vec3 baseColor {1.0f, 1.0f, 1.0f}; // diffuse color
    float metallic {0.0f};              // 0=dielectric, 1=metal
    float roughness {0.5f};             // 0=smooth,1=rough (used for spec exponent)
    float ior {1.45f};                  // index of refraction for dielectrics
    float alpha {1.0f};                 // opacity
    // For now, image textures are ignored or TODO; path can be stored if needed
    std::string baseColorTexturePath;
    std::shared_ptr<Texture> baseColorTexture;
};

} // namespace rt


