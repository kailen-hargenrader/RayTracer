// camera_to_rays.cpp
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <algorithm>
#include <sstream>

// Bring in camera implementation (Camera, Vec3, read_camera_from_json)
#include "../../s2899221/Code/persp_camera.cpp"

static void print_usage() {
    std::cerr << "Usage: camera_to_rays <cameras.json> <num_rays> <output.txt> [camera_id]\n";
    std::cerr << "- cameras.json: path to JSON with CAMERA->PERSP block (e.g., s2899221/ASCII/cameras.json)\n";
    std::cerr << "- num_rays: positive integer\n";
    std::cerr << "- output.txt: destination text file path\n";
    std::cerr << "- camera_id: optional id string; defaults to first camera if omitted\n";
}

// Extracts the first key name inside a JSON-like object block.
static bool extract_first_key(const std::string& block, std::string& out_key) {
    bool in_string = false; bool esc = false; bool seen_quote = false; std::string key;
    for (size_t i = 0; i < block.size(); ++i) {
        char c = block[i];
        if (in_string) {
            if (esc) { esc = false; key.push_back(c); }
            else if (c == '\\') esc = true;
            else if (c == '"') { in_string = false; seen_quote = true; }
            else key.push_back(c);
            continue;
        }
        if (c == '"') {
            if (!seen_quote) { in_string = true; key.clear(); }
        } else if (seen_quote && c == ':') {
            out_key = key; return true;
        }
    }
    return false;
}

// Note: extract_object_block is provided by the included persp_camera.cpp

int main(int argc, char** argv) {
    if (argc < 4) {
        print_usage();
        return 1;
    }

    const std::string json_path = argv[1];
    const int num_rays = std::max(0, std::stoi(argv[2]));
    const std::string out_path = argv[3];

    std::string camera_id;
    if (argc >= 5) {
        camera_id = argv[4];
    } else {
        // Discover first camera id from the JSON: CAMERA -> PERSP -> first key
        std::ifstream in(json_path);
        if (!in) { std::cerr << "Failed to open JSON file: " << json_path << "\n"; return 1; }
        std::stringstream buf; buf << in.rdbuf(); const std::string content = buf.str();
        std::string camera_block; if (!extract_object_block(content, "CAMERA", camera_block)) { std::cerr << "No CAMERA block\n"; return 1; }
        std::string persp_block; if (!extract_object_block(camera_block, "PERSP", persp_block)) { std::cerr << "No PERSP block\n"; return 1; }
        if (!extract_first_key(persp_block, camera_id)) { std::cerr << "No camera id found in PERSP block\n"; return 1; }
    }

    Camera cam;
    if (!read_camera_from_json(json_path, camera_id, cam)) {
        std::cerr << "Failed to read camera id '" << camera_id << "' from JSON\n";
        return 1;
    }

    int width = 0, height = 0;
    cam.getImageResolution(width, height);
    if (width <= 0 || height <= 0) {
        std::cerr << "Invalid camera resolution\n"; return 1;
    }

    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "Failed to open output: " << out_path << "\n";
        return 1;
    }

    std::mt19937 rng(123456u);
    std::uniform_real_distribution<double> dist_x(0.0, static_cast<double>(width));
    std::uniform_real_distribution<double> dist_y(0.0, static_cast<double>(height));

    // Header
    out << "# camera_id=" << camera_id << " width=" << width << " height=" << height << "\n";
    out << "# columns: ox oy oz dx dy dz\n";

    for (int i = 0; i < num_rays; ++i) {
        const double px = dist_x(rng);
        const double py = dist_y(rng);
        Vec3 origin, direction;
        cam.pixelToRay(px, py, origin, direction);
        out << origin.x << ' ' << origin.y << ' ' << origin.z << ' '
            << direction.x << ' ' << direction.y << ' ' << direction.z << "\n";
    }

    return 0;
}


