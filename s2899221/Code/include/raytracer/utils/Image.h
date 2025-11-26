#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "raytracer/math/Vec3.h"

namespace rt {

// Simple 8-bit RGB image with PPM read/write (P6 preferred)
class Image {
public:
    Image() = default;
    Image(int w, int h);

    void resize(int w, int h);
    int width() const { return m_width; }
    int height() const { return m_height; }

    // Access in [0,1] linear space
    void setPixel(int x, int y, const Vec3& rgb);
    Vec3 getPixel(int x, int y) const;

    // I/O
    bool writePPM(const std::string& path, bool binary = true) const;
    bool readPPM(const std::string& path); // supports P6 and P3

private:
    int m_width {0};
    int m_height {0};
    std::vector<uint8_t> m_pixels; // RGBRGB...
};

} // namespace rt


