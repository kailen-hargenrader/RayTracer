#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <regex>
#include <cmath>

#include "../../s2899221/Code/mesh.h"
#include "../../s2899221/Code/utils.h"

// Minimal JSON helpers reuse util_json

struct CameraCfg {
    double locx, locy, locz;
    double dirx, diry, dirz;
    double focal, sensor_w, sensor_h;
    int resx, resy;
};

static bool read_first_camera(const std::string& content, CameraCfg& out) {
    std::string cam_block;
    if (!util_json::extract_object_block(content, "CAMERA", cam_block)) return false;
    std::string persp_block;
    if (!util_json::extract_object_block(persp_block = cam_block, "PERSP", persp_block)) return false;

    // Extract first id by regex on quoted number keys
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(persp_block.begin(), persp_block.end(), id_rx);
    if (begin == std::sregex_iterator()) return false;
    const size_t key_pos = static_cast<size_t>(begin->position());
    std::string tail = persp_block.substr(key_pos);
    std::string one_cam; if (!util_json::extract_object_block(tail, begin->str(1), one_cam)) return false;

    double x=0,y=0,z=0;
    if (!util_json::parse_vec3(one_cam, "location", x, y, z)) return false; out.locx = x; out.locy = y; out.locz = z;
    if (!util_json::parse_vec3(one_cam, "direction", x, y, z)) return false; out.dirx = x; out.diry = y; out.dirz = z;
    if (!util_json::parse_number(one_cam, "focal_length", out.focal)) return false;
    if (!util_json::parse_number(one_cam, "sensor_width", out.sensor_w)) return false;
    if (!util_json::parse_number(one_cam, "sensor_height", out.sensor_h)) return false;
    if (!util_json::parse_vec2i(one_cam, "film_resolution", out.resx, out.resy)) return false;
    return true;
}

static bool read_meshes(const std::string& content, std::vector<Cube>& cubes, std::vector<Plane>& planes) {
    std::string mesh_block; if (!util_json::extract_object_block(content, "MESH", mesh_block)) return false;
    std::string cube_block; if (util_json::extract_object_block(mesh_block, "Cube", cube_block)) {
        auto v = Cube::read_from_json(cube_block); cubes.insert(cubes.end(), v.begin(), v.end());
    }
    std::string plane_block; if (util_json::extract_object_block(mesh_block, "Plane", plane_block)) {
        auto v = Plane::read_from_json(plane_block); planes.insert(planes.end(), v.begin(), v.end());
    }
    return true;
}

static void camera_basis(const double dir[3], double right[3], double up[3], double forward[3]) {
    // forward
    const double len = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    forward[0] = dir[0] / (len > 0 ? len : 1.0);
    forward[1] = dir[1] / (len > 0 ? len : 1.0);
    forward[2] = dir[2] / (len > 0 ? len : 1.0);
    // world up
    double wup[3] = {0,0,1};
    // right = normalize(cross(forward, wup))
    right[0] = forward[1]*wup[2] - forward[2]*wup[1];
    right[1] = forward[2]*wup[0] - forward[0]*wup[2];
    right[2] = forward[0]*wup[1] - forward[1]*wup[0];
    double rlen = std::sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]);
    if (rlen < 1e-9) { wup[0]=0; wup[1]=1; wup[2]=0; right[0] = forward[1]*wup[2] - forward[2]*wup[1]; right[1] = forward[2]*wup[0] - forward[0]*wup[2]; right[2] = forward[0]*wup[1] - forward[1]*wup[0]; rlen = std::sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]); }
    right[0]/=rlen; right[1]/=rlen; right[2]/=rlen;
    // up = cross(right, forward)
    up[0] = right[1]*forward[2] - right[2]*forward[1];
    up[1] = right[2]*forward[0] - right[0]*forward[2];
    up[2] = right[0]*forward[1] - right[1]*forward[0];
}

static void pixel_to_ray(const CameraCfg& cam, double px, double py, RayVec3& origin, RayVec3& direction) {
    double right[3], up[3], forward[3];
    double dir[3] = { cam.dirx, cam.diry, cam.dirz };
    camera_basis(dir, right, up, forward);
    const double nx = px / static_cast<double>(cam.resx);
    const double ny = py / static_cast<double>(cam.resy);
    const double x_sensor = cam.sensor_w * (nx - 0.5);
    const double y_sensor = cam.sensor_h * (0.5 - ny);
    double dw[3] = {
        right[0]*x_sensor + up[0]*y_sensor + forward[0]*cam.focal,
        right[1]*x_sensor + up[1]*y_sensor + forward[1]*cam.focal,
        right[2]*x_sensor + up[2]*y_sensor + forward[2]*cam.focal,
    };
    const double dlen = std::sqrt(dw[0]*dw[0] + dw[1]*dw[1] + dw[2]*dw[2]);
    origin = { cam.locx, cam.locy, cam.locz };
    direction = { dw[0]/dlen, dw[1]/dlen, dw[2]/dlen };
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: test_intersect_fn <path-to-json> <num-rays> [out.txt] [max_length]\n";
        return 1;
    }
    const std::string json_path = argv[1];
    const int num_rays = std::max(1, std::atoi(argv[2]));
    const std::string out_path = (argc >= 4) ? argv[3] : std::string("intersections.txt");
    const double max_len = (argc >= 5) ? std::max(0.0, std::atof(argv[4])) : 8.0;

    std::ifstream in(json_path);
    if (!in) { std::cerr << "Failed to open JSON: " << json_path << "\n"; return 1; }
    std::stringstream buffer; buffer << in.rdbuf();
    const std::string content = buffer.str();

    CameraCfg cam{};
    if (!read_first_camera(content, cam)) { std::cerr << "Failed to read camera from JSON\n"; return 1; }

    std::vector<Cube> cubes; std::vector<Plane> planes;
    read_meshes(content, cubes, planes);

    std::ofstream out(out_path);
    if (!out) { std::cerr << "Failed to open output: " << out_path << "\n"; return 1; }

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> ux(0.0, static_cast<double>(cam.resx));
    std::uniform_real_distribution<double> uy(0.0, static_cast<double>(cam.resy));

    for (int i = 0; i < num_rays; ++i) {
        double px = ux(rng), py = uy(rng);
        RayVec3 o, d; pixel_to_ray(cam, px, py, o, d);
        Ray ray(o, d);

        bool any_hit = false;
        double best_dist = 1e300;
        Hit best_hit;

        // Intersect cubes
        for (const auto& c : cubes) {
            Hit h;
            if (c.intersect(ray, h) && h.hasDistanceAlongRay()) {
                const double dist = h.getDistanceAlongRay();
                if (dist < best_dist) { best_dist = dist; best_hit = h; any_hit = true; }
            }
        }
        // TODO: planes intersect when Plane::intersect is implemented

        if (any_hit && best_dist <= max_len) {
            const auto& ip = best_hit.getIntersectionPoint();
            // Output ray segment from origin to hit (length = best_dist)
            out << o.x << ' ' << o.y << ' ' << o.z << ' ' << (ip.x - o.x) << ' ' << (ip.y - o.y) << ' ' << (ip.z - o.z) << ' ' << best_dist << '\n';
            // Output reflection ray with remaining length (max_len - best_dist)
            const double remain = std::max(0.0, max_len - best_dist);
            if (remain > 0.0) {
                const auto& refl = best_hit.getReflectedDirection();
                out << ip.x << ' ' << ip.y << ' ' << ip.z << ' ' << refl.x << ' ' << refl.y << ' ' << refl.z << ' ' << remain << '\n';
            }
        } else {
            // No hit within max_len: output primary ray capped at max_len and no reflection
            out << o.x << ' ' << o.y << ' ' << o.z << ' ' << d.x << ' ' << d.y << ' ' << d.z << ' ' << max_len << '\n';
        }
    }

    std::cout << "Wrote rays to " << out_path << "\n";
    return 0;
}


