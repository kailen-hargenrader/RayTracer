#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>

#include "../../s2899221/Code/mesh.h"

// Minimal helper: extract object block or array block following a key
static bool extract_block(const std::string& src, const std::string& key, std::string& out) {
    const std::string quoted_key = "\"" + key + "\"";
    size_t key_pos = src.find(quoted_key);
    if (key_pos == std::string::npos) return false;
    size_t colon_pos = src.find(':', key_pos + quoted_key.size());
    if (colon_pos == std::string::npos) return false;
    size_t i = colon_pos + 1;
    while (i < src.size() && std::isspace(static_cast<unsigned char>(src[i]))) ++i;
    if (i >= src.size()) return false;
    char open_char = src[i];
    if (open_char != '{' && open_char != '[') return false;
    char close_char = (open_char == '{') ? '}' : ']';
    int depth = 0; bool in_string = false; bool esc = false; size_t start = i;
    for (size_t j = i; j < src.size(); ++j) {
        char c = src[j];
        if (in_string) {
            if (esc) esc = false; else if (c == '\\') esc = true; else if (c == '"') in_string = false;
            continue;
        }
        if (c == '"') { in_string = true; continue; }
        if (c == open_char) { if (depth == 0) start = j; ++depth; }
        else if (c == close_char) { --depth; if (depth == 0) { out = src.substr(start + 1, j - start - 1); return true; } }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: json_to_console <path-to-json>\n";
        return 1;
    }

    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "Failed to open JSON file: " << argv[1] << "\n";
        return 1;
    }
    std::stringstream buffer; buffer << in.rdbuf();
    const std::string content = buffer.str();

    std::string mesh_block; if (!extract_block(content, "MESH", mesh_block)) {
        std::cerr << "No MESH block found.\n";
        return 1;
    }

    // Cube
    std::string cube_block; if (extract_block(mesh_block, "Cube", cube_block)) {
        auto cubes = Cube::read_from_json(cube_block);
        for (const auto& c : cubes) c.write_to_console(std::cout);
    }

    // Cylinder
    std::string cyl_block; if (extract_block(mesh_block, "Cylinder", cyl_block)) {
        auto cylinders = Cylinder::read_from_json(cyl_block);
        for (const auto& c : cylinders) c.write_to_console(std::cout);
    }

    // Sphere
    std::string sph_block; if (extract_block(mesh_block, "Sphere", sph_block)) {
        auto spheres = Sphere::read_from_json(sph_block);
        for (const auto& s : spheres) s.write_to_console(std::cout);
    }

    // Plane
    std::string pln_block; if (extract_block(mesh_block, "Plane", pln_block)) {
        auto planes = Plane::read_from_json(pln_block);
        for (const auto& p : planes) p.write_to_console(std::cout);
    }

    return 0;
}


