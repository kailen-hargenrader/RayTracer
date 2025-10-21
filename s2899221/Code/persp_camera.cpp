#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <cmath>
#include "utils.h"

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    Vec3 operator-(const Vec3& v) const {
        return Vec3(x - v.x, y - v.y, z - v.z);
    }

    Vec3 operator+(const Vec3& v) const {
        return Vec3(x + v.x, y + v.y, z + v.z);
    }

    Vec3 operator*(double s) const {
        return Vec3(x * s, y * s, z * s);
    }

    double length() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    Vec3 normalized() const {
        double len = length();
        if (len == 0) return Vec3(0,0,0);
        return Vec3(x/len, y/len, z/len);
    }
};

// Extracts the inner contents of an object block that follows a key in JSON-like text.
// Example: given src containing "\"FOO\": { ... }" and key "FOO",
// this returns the substring between the braces of that object (without the braces).
// JSON-like parsing helpers are implemented in utils.{h,cpp}

class Camera {
private:
    Vec3 position;
    double rotation[3][3];  // columns: right, up, forward
    double focal_length_mm;
    double sensor_width_mm;
    double sensor_height_mm;
    int image_width_px;
    int image_height_px;

    // Cross product implemented from scratch
    Vec3 cross(const Vec3& a, const Vec3& b) {
        return Vec3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }

public:
    Camera()
    : position(0,0,0), focal_length_mm(35.0), sensor_width_mm(36.0), sensor_height_mm(24.0),
      image_width_px(1920), image_height_px(1080)
    {
        // Initialize rotation matrix to identity
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                rotation[i][j] = (i == j) ? 1.0 : 0.0;
    }

    void setPosition(const Vec3& pos) {
        position = pos;
    }

    void setFocalLength(double focal_mm) {
        focal_length_mm = focal_mm;
    }

    void setSensorSize(double width_mm, double height_mm) {
        sensor_width_mm = width_mm;
        sensor_height_mm = height_mm;
    }

    void setImageResolution(int width_px, int height_px) {
        image_width_px = width_px;
        image_height_px = height_px;
    }

    // Expose the image resolution for sampling purposes.
    void getImageResolution(int& width_px, int& height_px) const {
        width_px = image_width_px;
        height_px = image_height_px;
    }

    // Update rotation matrix based on location and direction
    void updateRotationMatrix(const Vec3& cam_pos, const Vec3& cam_dir) {
        position = cam_pos;

        Vec3 forward = cam_dir.normalized();

        // Assume world up vector is (0,0,1)
        Vec3 world_up(0,0,1);

        // Use a right-handed basis: right = normalize(cross(forward, world_up))
        Vec3 right = cross(forward, world_up).normalized();

        // if right length is zero (camera looking straight up/down), fix by choosing another up
        if (right.length() < 1e-6) {
            // use (0,1,0) as up instead
            world_up = Vec3(0,1,0);
            right = cross(forward, world_up).normalized();
        }

        // up from right-handed basis
        Vec3 up = cross(right, forward);

        // Set rotation matrix columns: right, up, forward
        for (int i=0; i<3; ++i) {
            rotation[i][0] = (&right.x)[i];   // right vector
            rotation[i][1] = (&up.x)[i];      // up vector
            rotation[i][2] = (&forward.x)[i]; // forward vector
        }
    }

    // Project world point (x,y,z) to pixel (x,y)
    // Returns false if point behind camera, true otherwise
    bool worldToPixel(const Vec3& point_world, int& pixel_x, int& pixel_y) {
        // vector from camera to point
        Vec3 vec = point_world - position;

        // Transform vector to camera space: cam_coords = R^T * vec
        double Xc = rotation[0][0]*vec.x + rotation[1][0]*vec.y + rotation[2][0]*vec.z;
        double Yc = rotation[0][1]*vec.x + rotation[1][1]*vec.y + rotation[2][1]*vec.z;
        double Zc = rotation[0][2]*vec.x + rotation[1][2]*vec.y + rotation[2][2]*vec.z;

        if (Zc <= 0)
            return false; // point behind camera

        // Perspective projection on sensor plane (in mm)
        double x_sensor = focal_length_mm * (Xc / Zc);
        double y_sensor = focal_length_mm * (Yc / Zc);

        // Convert sensor coords to pixel coords
        pixel_x = static_cast<int>(((x_sensor / sensor_width_mm) + 0.5) * image_width_px);
        pixel_y = static_cast<int>((1.0 - ((y_sensor / sensor_height_mm) + 0.5)) * image_height_px);

        return true;
    }

    // Convert a pixel position (in pixel coordinates) to a world-space ray.
    // The ray originates at the camera position and points through the pixel center.
    void pixelToRay(double pixel_x, double pixel_y, Vec3& ray_origin, Vec3& ray_direction) const {
        // Normalize pixel coordinates to [0,1)
        const double nx = pixel_x / static_cast<double>(image_width_px);
        const double ny = pixel_y / static_cast<double>(image_height_px);

        // Map to sensor plane in millimeters
        const double x_sensor = sensor_width_mm * (nx - 0.5);
        const double y_sensor = sensor_height_mm * (0.5 - ny);

        // Direction in camera space
        const Vec3 dir_cam(x_sensor, y_sensor, focal_length_mm);

        // Rotate to world space: dir_world = R * dir_cam
        Vec3 dir_world;
        dir_world.x = rotation[0][0]*dir_cam.x + rotation[0][1]*dir_cam.y + rotation[0][2]*dir_cam.z;
        dir_world.y = rotation[1][0]*dir_cam.x + rotation[1][1]*dir_cam.y + rotation[1][2]*dir_cam.z;
        dir_world.z = rotation[2][0]*dir_cam.x + rotation[2][1]*dir_cam.y + rotation[2][2]*dir_cam.z;

        ray_origin = position;
        ray_direction = dir_world.normalized();
    }
};

// Reads camera parameters from a JSON file and updates the provided Camera instance.
// Returns true if a camera with the given id was found and applied; false otherwise.
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
