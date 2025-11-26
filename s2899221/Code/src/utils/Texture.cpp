#include "raytracer/utils/Texture.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>

namespace rt {

static std::vector<unsigned char> readAll(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    auto s = ss.str();
    return std::vector<unsigned char>(s.begin(), s.end());
}

static inline uint16_t readLE16(std::ifstream& in) {
    uint8_t b0 = 0, b1 = 0;
    in.read(reinterpret_cast<char*>(&b0), 1);
    in.read(reinterpret_cast<char*>(&b1), 1);
    return static_cast<uint16_t>(b0 | (static_cast<uint16_t>(b1) << 8));
}

static inline uint32_t readLE32(std::ifstream& in) {
    uint8_t b0 = 0, b1 = 0, b2 = 0, b3 = 0;
    in.read(reinterpret_cast<char*>(&b0), 1);
    in.read(reinterpret_cast<char*>(&b1), 1);
    in.read(reinterpret_cast<char*>(&b2), 1);
    in.read(reinterpret_cast<char*>(&b3), 1);
    return static_cast<uint32_t>(b0)
        | (static_cast<uint32_t>(b1) << 8)
        | (static_cast<uint32_t>(b2) << 16)
        | (static_cast<uint32_t>(b3) << 24);
}

static bool loadBMP24or32(const std::string& path, int& outW, int& outH, std::vector<unsigned char>& outRGB) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    // BITMAPFILEHEADER (14 bytes)
    uint8_t magic[2] = {0,0};
    in.read(reinterpret_cast<char*>(magic), 2);
    if (magic[0] != 'B' || magic[1] != 'M') return false;
    uint32_t fileSize = readLE32(in);
    (void)fileSize;
    uint16_t res1 = readLE16(in);
    uint16_t res2 = readLE16(in);
    (void)res1; (void)res2;
    uint32_t pixelOffset = readLE32(in);

    // DIB header - assume BITMAPINFOHEADER (40 bytes)
    uint32_t dibSize = readLE32(in);
    if (dibSize < 40) return false;
    int32_t width = static_cast<int32_t>(readLE32(in));
    int32_t height = static_cast<int32_t>(readLE32(in));
    uint16_t planes = readLE16(in);
    uint16_t bitCount = readLE16(in);
    uint32_t compression = readLE32(in);
    uint32_t sizeImage = readLE32(in);
    uint32_t xppm = readLE32(in);
    uint32_t yppm = readLE32(in);
    uint32_t clrUsed = readLE32(in);
    uint32_t clrImportant = readLE32(in);
    (void)sizeImage; (void)xppm; (void)yppm; (void)clrUsed; (void)clrImportant;

    if (planes != 1) return false;
    if (compression != 0) return false; // BI_RGB only
    if (!(bitCount == 24 || bitCount == 32)) return false;
    if (width <= 0 || height == 0) return false;

    const bool topDown = (height < 0);
    const int h = std::abs(height);
    const int w = width;
    const int bytesPerPixel = bitCount / 8;
    const int rowStrideFile = ((w * bytesPerPixel + 3) / 4) * 4; // padded to multiple of 4

    // Seek to pixel data
    in.seekg(static_cast<std::streamoff>(pixelOffset), std::ios::beg);
    if (!in) return false;

    std::vector<unsigned char> rgb;
    rgb.resize(static_cast<size_t>(w) * static_cast<size_t>(h) * 3u);
    std::vector<uint8_t> rowBuf(static_cast<size_t>(rowStrideFile));

    for (int row = 0; row < h; ++row) {
        int dstRow = topDown ? row : (h - 1 - row);
        in.read(reinterpret_cast<char*>(rowBuf.data()), rowStrideFile);
        if (in.gcount() != rowStrideFile) return false;
        for (int x = 0; x < w; ++x) {
            size_t srcIdx = static_cast<size_t>(x) * static_cast<size_t>(bytesPerPixel);
            size_t dstIdx = (static_cast<size_t>(dstRow) * static_cast<size_t>(w) + static_cast<size_t>(x)) * 3u;
            uint8_t b = rowBuf[srcIdx + 0];
            uint8_t g = rowBuf[srcIdx + 1];
            uint8_t r = rowBuf[srcIdx + 2];
            rgb[dstIdx + 0] = r;
            rgb[dstIdx + 1] = g;
            rgb[dstIdx + 2] = b;
        }
    }

    outW = w;
    outH = h;
    outRGB.swap(rgb);
    return true;
}

bool Texture::loadFromFile(const std::string& path) {
    m_width = 0;
    m_height = 0;
    m_channels = 0;
    m_pixels.clear();

    std::filesystem::path p(path);
    std::string ext = p.has_extension() ? p.extension().string() : std::string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    int w = 0, h = 0;
    std::vector<unsigned char> rgb;

    // Support BMP (24/32-bit, uncompressed). This requires only the C++ standard library.
    if (ext == ".bmp") {
        if (!loadBMP24or32(path, w, h, rgb)) return false;
    } else {
        // Unknown extension: do not attempt to decode formats that require external libraries.
        // Return false to allow caller to try alternative paths (e.g., .bmp or .ppm).
        return false;
    }

    m_width = w;
    m_height = h;
    m_channels = 3;
    m_pixels.swap(rgb);
    return true;
}

Vec3 Texture::sample(float u, float v) const {
    if (!isValid()) return {1,1,1};
    // Wrap repeat
    u = u - std::floor(u);
    v = v - std::floor(v);
    int ix = std::max(0, std::min(m_width - 1, static_cast<int>(u * m_width)));
    int iy = std::max(0, std::min(m_height - 1, static_cast<int>((1.0f - v) * m_height)));
    size_t idx = (static_cast<size_t>(iy) * static_cast<size_t>(m_width) + static_cast<size_t>(ix)) * static_cast<size_t>(m_channels);
	float r_srgb = static_cast<float>(m_pixels[idx + 0]) / 255.0f;
	float g_srgb = static_cast<float>(m_pixels[idx + 1]) / 255.0f;
	float b_srgb = static_cast<float>(m_pixels[idx + 2]) / 255.0f;

	// Convert from sRGB to linear for correct lighting computations
	auto srgbToLinear = [](float c) -> float {
		return (c <= 0.04045f) ? (c / 12.92f)
			: std::pow((c + 0.055f) / 1.055f, 2.4f);
	};
	float r = srgbToLinear(r_srgb);
	float g = srgbToLinear(g_srgb);
	float b = srgbToLinear(b_srgb);
	return {r, g, b};
}

} // namespace rt


