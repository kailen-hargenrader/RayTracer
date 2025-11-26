#pragma once

#include <string>
#include <vector>
#include "raytracer/math/Vec3.h"

namespace rt {

class Texture {
public:
    Texture() = default;
    ~Texture() = default;

    bool loadFromFile(const std::string& path);
    bool isValid() const { return m_width > 0 && m_height > 0 && !m_pixels.empty(); }
    Vec3 sample(float u, float v) const; // nearest for simplicity

private:
    int m_width {0};
    int m_height {0};
    int m_channels {0}; // 3 or 4
    std::vector<unsigned char> m_pixels; // row-major, RGB/RGBA
};

} // namespace rt


