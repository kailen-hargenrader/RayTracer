#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <cmath>

#include "persp_camera.h"
#include "utils.h"

Vec3::Vec3() : x(0), y(0), z(0) {}
Vec3::Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

Vec3 Vec3::operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
Vec3 Vec3::operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
Vec3 Vec3::operator*(double s) const { return Vec3(x * s, y * s, z * s); }
double Vec3::length() const { return std::sqrt(x*x + y*y + z*z); }
Vec3 Vec3::normalized() const { double len = length(); if (len == 0) return Vec3(0,0,0); return Vec3(x/len, y/len, z/len); }

Camera::Camera()
    : position(0,0,0), focal_length_mm(35.0), sensor_width_mm(36.0), sensor_height_mm(24.0),
      image_width_px(1920), image_height_px(1080) {
    for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) rotation[i][j] = (i == j) ? 1.0 : 0.0;
}

void Camera::setPosition(const Vec3& pos) { position = pos; }
void Camera::setFocalLength(double focal_mm) { focal_length_mm = focal_mm; }
void Camera::setSensorSize(double width_mm, double height_mm) { sensor_width_mm = width_mm; sensor_height_mm = height_mm; }
void Camera::setImageResolution(int width_px, int height_px) { image_width_px = width_px; image_height_px = height_px; }
void Camera::getImageResolution(int& width_px, int& height_px) const { width_px = image_width_px; height_px = image_height_px; }

Vec3 Camera::cross(const Vec3& a, const Vec3& b) {
    return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

void Camera::updateRotationMatrix(const Vec3& cam_pos, const Vec3& cam_dir) {
    position = cam_pos;
    Vec3 forward = cam_dir.normalized();
    Vec3 world_up(0,0,1);
    Vec3 right = cross(forward, world_up).normalized();
    if (right.length() < 1e-6) { world_up = Vec3(0,1,0); right = cross(forward, world_up).normalized(); }
    Vec3 up = cross(right, forward);
    for (int i=0; i<3; ++i) {
        rotation[i][0] = (&right.x)[i];
        rotation[i][1] = (&up.x)[i];
        rotation[i][2] = (&forward.x)[i];
    }
}

void Camera::updateRotationMatrix(const Vec3& cam_pos, const Vec3& cam_dir, const Vec3& cam_up) {
    position = cam_pos;
    Vec3 forward = cam_dir.normalized();
    // Use provided up to preserve roll, but re-orthonormalize
    Vec3 up_guess = cam_up.normalized();
    // If up is degenerate or parallel to forward, fall back to world-up
    if (up_guess.length() < 1e-6 || std::abs(forward.x*up_guess.x + forward.y*up_guess.y + forward.z*up_guess.z) > 0.999) {
        updateRotationMatrix(cam_pos, cam_dir);
        return;
    }
    Vec3 right = cross(forward, up_guess).normalized();
    Vec3 up = cross(right, forward);
    for (int i=0; i<3; ++i) {
        rotation[i][0] = (&right.x)[i];
        rotation[i][1] = (&up.x)[i];
        rotation[i][2] = (&forward.x)[i];
    }
}
bool Camera::worldToPixel(const Vec3& point_world, int& pixel_x, int& pixel_y) {
    Vec3 vec = point_world - position;
    double Xc = rotation[0][0]*vec.x + rotation[1][0]*vec.y + rotation[2][0]*vec.z;
    double Yc = rotation[0][1]*vec.x + rotation[1][1]*vec.y + rotation[2][1]*vec.z;
    double Zc = rotation[0][2]*vec.x + rotation[1][2]*vec.y + rotation[2][2]*vec.z;
    if (Zc <= 0) return false;
    // Match Blender's sensor fit behavior: fit horizontally if image aspect >= sensor aspect, else fit vertically
    const double aspect_img = static_cast<double>(image_width_px) / static_cast<double>(image_height_px);
    const double aspect_sensor = (sensor_height_mm > 0.0) ? (sensor_width_mm / sensor_height_mm) : aspect_img;
    double eff_w = sensor_width_mm;
    double eff_h = sensor_height_mm;
    if (aspect_img >= aspect_sensor) {
        // Horizontal fit: width fixed, height derived from image aspect
        eff_w = sensor_width_mm;
        eff_h = eff_w / aspect_img;
    } else {
        // Vertical fit: height fixed, width derived
        eff_h = sensor_height_mm;
        eff_w = eff_h * aspect_img;
    }
    double x_sensor = focal_length_mm * (Xc / Zc);
    double y_sensor = focal_length_mm * (Yc / Zc);
    pixel_x = static_cast<int>(((x_sensor / eff_w) + 0.5) * image_width_px);
    pixel_y = static_cast<int>((1.0 - ((y_sensor / eff_h) + 0.5)) * image_height_px);
    return true;
}

void Camera::pixelToRay(double pixel_x, double pixel_y, Vec3& ray_origin, Vec3& ray_direction) const {
    const double nx = pixel_x / static_cast<double>(image_width_px);
    const double ny = pixel_y / static_cast<double>(image_height_px);
    // Effective sensor size to match Blender's sensor fit rules
    const double aspect_img = static_cast<double>(image_width_px) / static_cast<double>(image_height_px);
    const double aspect_sensor = (sensor_height_mm > 0.0) ? (sensor_width_mm / sensor_height_mm) : aspect_img;
    double eff_w = sensor_width_mm;
    double eff_h = sensor_height_mm;
    if (aspect_img >= aspect_sensor) {
        eff_w = sensor_width_mm;
        eff_h = eff_w / aspect_img;
    } else {
        eff_h = sensor_height_mm;
        eff_w = eff_h * aspect_img;
    }
    const double x_sensor = eff_w * (nx - 0.5);
    const double y_sensor = eff_h * (0.5 - ny);
    const Vec3 dir_cam(x_sensor, y_sensor, focal_length_mm);
    Vec3 dir_world;
    dir_world.x = rotation[0][0]*dir_cam.x + rotation[0][1]*dir_cam.y + rotation[0][2]*dir_cam.z;
    dir_world.y = rotation[1][0]*dir_cam.x + rotation[1][1]*dir_cam.y + rotation[1][2]*dir_cam.z;
    dir_world.z = rotation[2][0]*dir_cam.x + rotation[2][1]*dir_cam.y + rotation[2][2]*dir_cam.z;
    ray_origin = position;
    ray_direction = dir_world.normalized();
}

bool read_camera_from_json(const std::string& json_filepath, const std::string& camera_id, Camera& camera) {
    std::ifstream in(json_filepath);
    if (!in) {
        std::cerr << "Failed to open JSON file: " << json_filepath << "\n";
        return false;
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string content = buffer.str();

    // Extract blocks step-by-step: CAMERA -> PERSP -> camera_id
    std::string camera_block;
    if (!util_json::extract_object_block(content, "CAMERA", camera_block)) {
        return false;
    }
    std::string persp_block;
    if (!util_json::extract_object_block(camera_block, "PERSP", persp_block)) {
        return false;
    }
    std::string cam_block;
    if (!util_json::extract_object_block(persp_block, camera_id, cam_block)) {
        return false;
    }

    auto extract_number = [](const std::string& src, const std::string& key, double& out) -> bool {
        std::regex r("\\\"" + key + "\\\"\\s*:\\s*([-+eE0-9\\.]+)");
        std::smatch m; if (!std::regex_search(src, m, r)) return false; out = std::stod(m[1]); return true; };

    Vec3 location, direction;
    double focal = 0.0, sensor_w = 0.0, sensor_h = 0.0; int resx = 0, resy = 0;
    bool ok = true;
    {
        double x=0,y=0,z=0; ok &= util_json::parse_vec3(cam_block, "location", x, y, z); location = Vec3(x,y,z);
    }
    {
        double x=0,y=0,z=0; ok &= util_json::parse_vec3(cam_block, "direction", x, y, z); direction = Vec3(x,y,z);
    }
    ok &= util_json::parse_number(cam_block, "focal_length", focal);
    ok &= util_json::parse_number(cam_block, "sensor_width", sensor_w);
    ok &= util_json::parse_number(cam_block, "sensor_height", sensor_h);
    ok &= util_json::parse_vec2i(cam_block, "film_resolution", resx, resy);
    if (!ok) return false;

    camera.setFocalLength(focal);
    camera.setSensorSize(sensor_w, sensor_h);
    camera.setImageResolution(resx, resy);
    camera.updateRotationMatrix(location, direction);
    return true;
}
