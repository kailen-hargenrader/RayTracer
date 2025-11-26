#include "raytracer/utils/Texture.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

namespace rt {

static std::vector<unsigned char> readAll(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    auto s = ss.str();
    return std::vector<unsigned char>(s.begin(), s.end());
}

bool Texture::loadFromFile(const std::string& path) {
    // Image loader disabled (no external deps). Return false to fallback to base color.
    (void)path;
    return false;
}

Vec3 Texture::sample(float u, float v) const {
    if (!isValid()) return {1,1,1};
    // Wrap repeat
    u = u - std::floor(u);
    v = v - std::floor(v);
    int ix = std::max(0, std::min(m_width - 1, static_cast<int>(u * m_width)));
    int iy = std::max(0, std::min(m_height - 1, static_cast<int>((1.0f - v) * m_height)));
    size_t idx = (static_cast<size_t>(iy) * static_cast<size_t>(m_width) + static_cast<size_t>(ix)) * static_cast<size_t>(m_channels);
    float r = static_cast<float>(m_pixels[idx + 0]) / 255.0f;
    float g = static_cast<float>(m_pixels[idx + 1]) / 255.0f;
    float b = static_cast<float>(m_pixels[idx + 2]) / 255.0f;
    return {r, g, b};
}

} // namespace rt


