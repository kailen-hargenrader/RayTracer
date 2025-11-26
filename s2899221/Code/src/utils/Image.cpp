#include "raytracer/utils/Image.h"
#include <fstream>
#include <sstream>
#include <algorithm>

namespace rt {

Image::Image(int w, int h) {
    resize(w, h);
}

void Image::resize(int w, int h) {
    m_width = std::max(0, w);
    m_height = std::max(0, h);
    m_pixels.assign(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3u, 0);
}

void Image::setPixel(int x, int y, const Vec3& rgb) {
    if (x < 0 || y < 0 || x >= m_width || y >= m_height) return;
	// Convert from linear to sRGB for perceptual display/output
	auto linearToSrgb = [](float c) -> float {
		c = std::max(0.0f, std::min(1.0f, c));
		return (c <= 0.0031308f) ? (12.92f * c)
			: (1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f);
	};
	float r_srgb = linearToSrgb(rgb.x);
	float g_srgb = linearToSrgb(rgb.y);
	float b_srgb = linearToSrgb(rgb.z);
	uint8_t r = static_cast<uint8_t>(std::round(r_srgb * 255.0f));
	uint8_t g = static_cast<uint8_t>(std::round(g_srgb * 255.0f));
	uint8_t b = static_cast<uint8_t>(std::round(b_srgb * 255.0f));
    size_t idx = (static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)) * 3u;
    m_pixels[idx + 0] = r;
    m_pixels[idx + 1] = g;
    m_pixels[idx + 2] = b;
}

Vec3 Image::getPixel(int x, int y) const {
    if (x < 0 || y < 0 || x >= m_width || y >= m_height) return {0,0,0};
    size_t idx = (static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)) * 3u;
    return {
        static_cast<float>(m_pixels[idx + 0]) / 255.0f,
        static_cast<float>(m_pixels[idx + 1]) / 255.0f,
        static_cast<float>(m_pixels[idx + 2]) / 255.0f
    };
}

bool Image::writePPM(const std::string& path, bool binary) const {
    if (m_width <= 0 || m_height <= 0) return false;
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;

    if (binary) {
        out << "P6\n" << m_width << " " << m_height << "\n255\n";
        out.write(reinterpret_cast<const char*>(m_pixels.data()), static_cast<std::streamsize>(m_pixels.size()));
    } else {
        out << "P3\n" << m_width << " " << m_height << "\n255\n";
        for (int y = 0; y < m_height; ++y) {
            for (int x = 0; x < m_width; ++x) {
                size_t idx = (static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)) * 3u;
                out << static_cast<int>(m_pixels[idx + 0]) << " "
                    << static_cast<int>(m_pixels[idx + 1]) << " "
                    << static_cast<int>(m_pixels[idx + 2]) << "\n";
            }
        }
    }
    return true;
}

static bool readToken(std::istream& in, std::string& token) {
    token.clear();
    while (in >> token) {
        if (!token.empty() && token[0] == '#') {
            std::string rest;
            std::getline(in, rest);
            continue;
        }
        return true;
    }
    return false;
}

bool Image::readPPM(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    std::string magic;
    in >> magic;
    if (magic != "P6" && magic != "P3") return false;

    std::string tok;
    if (!readToken(in, tok)) return false;
    int w = std::stoi(tok);
    if (!readToken(in, tok)) return false;
    int h = std::stoi(tok);
    if (!readToken(in, tok)) return false;
    int maxv = std::stoi(tok);
    if (maxv <= 0 || maxv > 65535) return false;
    resize(w, h);

    if (magic == "P6") {
        // Consume single whitespace after header
        in.get();
        size_t expected = static_cast<size_t>(w) * static_cast<size_t>(h) * 3u;
        m_pixels.assign(expected, 0);
        in.read(reinterpret_cast<char*>(m_pixels.data()), static_cast<std::streamsize>(expected));
        return in.good() || in.eof();
    } else {
        // ASCII
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                std::string r,g,b;
                if (!readToken(in, r) || !readToken(in, g) || !readToken(in, b)) return false;
                int ri = std::stoi(r);
                int gi = std::stoi(g);
                int bi = std::stoi(b);
                size_t idx = (static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)) * 3u;
                auto clamp8 = [](int v){ return static_cast<uint8_t>(std::max(0, std::min(255, v))); };
                if (maxv != 255) {
                    // Normalize to 255
                    ri = static_cast<int>(std::round((ri / static_cast<float>(maxv)) * 255.0f));
                    gi = static_cast<int>(std::round((gi / static_cast<float>(maxv)) * 255.0f));
                    bi = static_cast<int>(std::round((bi / static_cast<float>(maxv)) * 255.0f));
                }
                m_pixels[idx+0] = clamp8(ri);
                m_pixels[idx+1] = clamp8(gi);
                m_pixels[idx+2] = clamp8(bi);
            }
        }
        return true;
    }
}

} // namespace rt


