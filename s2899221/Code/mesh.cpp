#include "mesh.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <cctype>

// Utilities copied/adapted from persp_camera.cpp for simple JSON-like parsing
static bool extract_object_block(const std::string& src, const std::string& key, std::string& out_block) {
    const std::string quoted_key = "\"" + key + "\"";
    const size_t key_pos = src.find(quoted_key);
    if (key_pos == std::string::npos) return false;

    size_t colon_pos = src.find(':', key_pos + quoted_key.size());
    if (colon_pos == std::string::npos) return false;

    // Allow either object { } or array [ ] blocks
    size_t opener_pos = std::string::npos;
    char open_char = 0;
    char close_char = 0;
    // Find next non-space char after colon
    size_t i = colon_pos + 1;
    while (i < src.size() && std::isspace(static_cast<unsigned char>(src[i]))) ++i;
    if (i >= src.size()) return false;
    if (src[i] == '{' || src[i] == '[') {
        opener_pos = i;
        open_char = src[i];
        close_char = (open_char == '{') ? '}' : ']';
    } else {
        // Not a block
        return false;
    }

    int depth = 0;
    bool in_string = false;
    bool is_escaped = false;
    for (size_t j = opener_pos; j < src.size(); ++j) {
        const char c = src[j];
        if (in_string) {
            if (is_escaped) {
                is_escaped = false;
            } else if (c == '\\') {
                is_escaped = true;
            } else if (c == '"') {
                in_string = false;
            }
            continue;
        }
        if (c == '"') { in_string = true; continue; }
        if (c == open_char) {
            if (depth == 0) opener_pos = j;
            ++depth;
        } else if (c == close_char) {
            --depth;
            if (depth == 0) {
                out_block = src.substr(opener_pos + 1, j - opener_pos - 1);
                return true;
            }
        }
    }
    return false;
}

static bool parse_number(const std::string& src, const std::string& key, double& out) {
    std::regex r("\\\"" + key + "\\\"\\s*:\\s*([-+eE0-9\\.]+)");
    std::smatch m; if (!std::regex_search(src, m, r)) return false; out = std::stod(m[1]); return true;
}

static bool parse_vec3(const std::string& src, const std::string& key, Float3& out) {
    std::string block; if (!extract_object_block(src, key, block)) return false;
    double x, y, z; if (!parse_number(block, "x", x)) return false; if (!parse_number(block, "y", y)) return false; if (!parse_number(block, "z", z)) return false;
    out = Float3{ x, y, z }; return true;
}

// ---------- Cube ----------
Cube::Cube() : m_translation{0,0,0}, m_rotation{0,0,0}, m_scale(1.0) {}
Cube::Cube(const Float3& translation, const EulerAngles& rotation, double scale)
    : m_translation(translation), m_rotation(rotation), m_scale(scale) {}

std::vector<Cube> Cube::read_from_json(const std::string& class_block) {
    // class_block expected to contain entries keyed by id: { "1": { ... }, "2": { ... }, ... }
    std::vector<Cube> result;

    // Find each object id by regex on quoted number keys
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        // Extract the object block following this id
        std::string sub;
        // Build a temporary string starting from this position to ensure extract works
        std::string tail = class_block.substr(key_pos);
        if (!extract_object_block(tail, it->str(1), sub)) continue;

        Float3 tr{0,0,0}; EulerAngles rot{0,0,0}; double scale = 1.0;
        parse_vec3(sub, "translation", tr);
        // rotation fields named roll, pitch, yaw
        double r=0,p=0,y=0; parse_number(sub, "roll", r); parse_number(sub, "pitch", p); parse_number(sub, "yaw", y);
        rot = EulerAngles{ r, p, y };
        parse_number(sub, "scale", scale);

        result.emplace_back(tr, rot, scale);
    }
    return result;
}

void Cube::write_to_console(std::ostream& out) const {
    out << "Cube(translation=[" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z
        << "], rotation(rpy)=[" << m_rotation.roll << ", " << m_rotation.pitch << ", " << m_rotation.yaw
        << "], scale=" << m_scale << ")\n";
}

// ---------- Cylinder ----------
Cylinder::Cylinder() : m_translation{0,0,0}, m_rotation{0,0,0}, m_scale(1.0), m_radius(1.0), m_length(1.0) {}
Cylinder::Cylinder(const Float3& translation, const EulerAngles& rotation, double scale, double radius, double length)
    : m_translation(translation), m_rotation(rotation), m_scale(scale), m_radius(radius), m_length(length) {}

std::vector<Cylinder> Cylinder::read_from_json(const std::string& class_block) {
    std::vector<Cylinder> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!extract_object_block(tail, it->str(1), sub)) continue;

        Float3 tr{0,0,0}; EulerAngles rot{0,0,0}; double scale = 1.0; double radius = 1.0; double length = 1.0;
        parse_vec3(sub, "translation", tr);
        double r=0,p=0,y=0; parse_number(sub, "roll", r); parse_number(sub, "pitch", p); parse_number(sub, "yaw", y);
        rot = EulerAngles{ r, p, y };
        parse_number(sub, "scale", scale);
        parse_number(sub, "radius", radius);
        parse_number(sub, "length", length);

        result.emplace_back(tr, rot, scale, radius, length);
    }
    return result;
}

void Cylinder::write_to_console(std::ostream& out) const {
    out << "Cylinder(translation=[" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z
        << "], rotation(rpy)=[" << m_rotation.roll << ", " << m_rotation.pitch << ", " << m_rotation.yaw
        << "], scale=" << m_scale << ", radius=" << m_radius << ", length=" << m_length << ")\n";
}

// ---------- Sphere ----------
Sphere::Sphere() : m_location{0,0,0}, m_radius(1.0) {}
Sphere::Sphere(const Float3& location, double radius) : m_location(location), m_radius(radius) {}

std::vector<Sphere> Sphere::read_from_json(const std::string& class_block) {
    std::vector<Sphere> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!extract_object_block(tail, it->str(1), sub)) continue;

        Float3 loc{0,0,0}; double radius = 1.0;
        parse_vec3(sub, "location", loc);
        parse_number(sub, "radius", radius);
        result.emplace_back(loc, radius);
    }
    return result;
}

void Sphere::write_to_console(std::ostream& out) const {
    out << "Sphere(location=[" << m_location.x << ", " << m_location.y << ", " << m_location.z
        << "], radius=" << m_radius << ")\n";
}

// ---------- Plane ----------
Plane::Plane() : m_corners() {}
Plane::Plane(const std::vector<Float3>& corners) : m_corners(corners) {}

std::vector<Plane> Plane::read_from_json(const std::string& class_block) {
    std::vector<Plane> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!extract_object_block(tail, it->str(1), sub)) continue;

        // corners is an array of 4 objects {x,y,z}
        std::string corners_block;
        if (!extract_object_block(sub, "corners", corners_block)) {
            continue;
        }
        // corners_block here is content between '[' and ']' thanks to support in extract_object_block
        // Parse up to 4 corner objects in order of appearance
        std::vector<Float3> corners;
        size_t search_offset = 0;
        while (true) {
            const size_t brace_pos = corners_block.find('{', search_offset);
            if (brace_pos == std::string::npos) break;

            // Manual brace match to find end of this corner object
            int depth = 0; bool in_str = false; bool esc = false; size_t end_pos = std::string::npos;
            for (size_t j = brace_pos; j < corners_block.size(); ++j) {
                char c = corners_block[j];
                if (in_str) { if (esc) esc = false; else if (c == '\\') esc = true; else if (c == '"') in_str = false; continue; }
                if (c == '"') { in_str = true; continue; }
                if (c == '{') { ++depth; }
                else if (c == '}') { --depth; if (depth == 0) { end_pos = j; break; } }
            }
            if (end_pos == std::string::npos) break;

            const std::string one_corner = corners_block.substr(brace_pos + 1, end_pos - brace_pos - 1);
            search_offset = end_pos + 1;

            double x=0,y=0,z=0; parse_number(one_corner, "x", x); parse_number(one_corner, "y", y); parse_number(one_corner, "z", z);
            corners.push_back(Float3{ x, y, z });
            if (corners.size() == 4) break;
        }
        if (!corners.empty()) {
            result.emplace_back(corners);
        }
    }
    return result;
}

void Plane::write_to_console(std::ostream& out) const {
    out << "Plane(corners=[";
    for (size_t i = 0; i < m_corners.size(); ++i) {
        const auto& c = m_corners[i];
        out << "[" << c.x << ", " << c.y << ", " << c.z << "]";
        if (i + 1 < m_corners.size()) out << ", ";
    }
    out << "])\n";
}


