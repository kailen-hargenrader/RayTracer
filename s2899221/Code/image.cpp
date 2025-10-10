#include "image.h"

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cctype>

static void skipComments(std::istream& in) {
    // Skip whitespace and comments starting with '#'
    int c = in.peek();
    while (std::isspace(c) || c == '#') {
        if (c == '#') {
            std::string dummy;
            std::getline(in, dummy);
        } else {
            in.get();
        }
        c = in.peek();
    }
}

Image::Image(const std::string& filename)
    : m_width(0), m_height(0), m_maxValue(255) {
    loadFromFile(filename);
}

Image::Image(int width, int height)
    : m_width(width), m_height(height), m_maxValue(255), m_data(static_cast<size_t>(width) * static_cast<size_t>(height)) {}

int Image::width() const { return m_width; }
int Image::height() const { return m_height; }
int Image::maxValue() const { return m_maxValue; }

Pixel Image::getPixel(int x, int y) const {
    return m_data[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)];
}

static inline unsigned char clamp_to_max(int v, int maxv) {
    if (v < 0) return 0;
    if (v > maxv) return static_cast<unsigned char>(maxv);
    return static_cast<unsigned char>(v);
}

void Image::setPixel(int x, int y, int r, int g, int b) {
    Pixel p;
    const int maxv = m_maxValue;
    p.r = clamp_to_max(r, maxv);
    p.g = clamp_to_max(g, maxv);
    p.b = clamp_to_max(b, maxv);
    m_data[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] = p;
}

void Image::addToPixel(int x, int y, int dr, int dg, int db) {
    Pixel p = getPixel(x, y);
    const int maxv = m_maxValue;
    p.r = clamp_to_max(static_cast<int>(p.r) + dr, maxv);
    p.g = clamp_to_max(static_cast<int>(p.g) + dg, maxv);
    p.b = clamp_to_max(static_cast<int>(p.b) + db, maxv);
    m_data[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] = p;
}

void Image::setMaxValue(int maxVal) {
    // Only 8-bit PPM is supported here; clamp to [1,255]
    if (maxVal < 1 || maxVal > 255) throw std::runtime_error("Invalid PPM max value, must be [1,255].");
    m_maxValue = maxVal;
}

void Image::write(const std::string& filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    // 1..255 supported; write using current max value
    out << "P6\n" << m_width << " " << m_height << "\n" << m_maxValue << "\n";
    // Data is already in 0..maxval; clamp once and write as-is
    std::vector<unsigned char> buf;
    buf.reserve(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3u);
    for (const Pixel& p : m_data) {
        buf.push_back(clamp_to_max(p.r, m_maxValue));
        buf.push_back(clamp_to_max(p.g, m_maxValue));
        buf.push_back(clamp_to_max(p.b, m_maxValue));
    }
    out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
    if (!out) {
        throw std::runtime_error("Failed while writing image data to: " + filename);
    }
}

void Image::loadFromFile(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open image file: " + filename);
    }

    // Read magic number
    std::string magic;
    in >> magic;
    if (magic != "P6" && magic != "P3") {
        throw std::runtime_error("Unsupported PPM format (expected P3 or P6): " + magic);
    }

    skipComments(in);
    in >> m_width; skipComments(in);
    in >> m_height; skipComments(in);
    int maxval = 0; in >> maxval; in.get(); // consume one whitespace/newline after maxval
    // Only 8-bit max values are supported (1..255). 16-bit PPM is not supported.
    if (m_width <= 0 || m_height <= 0 || maxval <= 0) {
        throw std::runtime_error("Invalid PPM header values");
    }
    if (maxval > 255) {
        throw std::runtime_error("Unsupported PPM max value (>255 not supported)");
    }

    m_maxValue = maxval;
    m_data.resize(static_cast<size_t>(m_width) * static_cast<size_t>(m_height));

    if (magic == "P6") {
        // Binary RGB packed immediately after header
        const size_t numBytes = static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3u;
        std::vector<unsigned char> buf(numBytes);
        in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
        if (static_cast<size_t>(in.gcount()) != numBytes) {
            throw std::runtime_error("Unexpected EOF while reading binary PPM data");
        }
        for (size_t i = 0, p = 0; i < m_data.size(); ++i) {
            m_data[i].r = buf[p++];
            m_data[i].g = buf[p++];
            m_data[i].b = buf[p++];
        }
    } else {
        // P3 ASCII
        for (size_t i = 0; i < m_data.size(); ++i) {
            int r, g, b;
            skipComments(in);
            if (!(in >> r)) throw std::runtime_error("Unexpected EOF reading P3 data (r)");
            skipComments(in);
            if (!(in >> g)) throw std::runtime_error("Unexpected EOF reading P3 data (g)");
            skipComments(in);
            if (!(in >> b)) throw std::runtime_error("Unexpected EOF reading P3 data (b)");
            if (r < 0 || r > maxval || g < 0 || g > maxval || b < 0 || b > maxval) {
                throw std::runtime_error("P3 pixel out of range");
            }
            m_data[i].r = static_cast<unsigned char>(r);
            m_data[i].g = static_cast<unsigned char>(g);
            m_data[i].b = static_cast<unsigned char>(b);
        }
    }
}


