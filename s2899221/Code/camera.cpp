#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <cmath>

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

    // Update rotation matrix based on location and direction
    void updateRotationMatrix(const Vec3& cam_pos, const Vec3& cam_dir) {
        position = cam_pos;

        Vec3 forward = cam_dir.normalized();

        // Assume world up vector is (0,0,1)
        Vec3 world_up(0,0,1);

        // right = normalize(cross(world_up, forward))
        Vec3 right = cross(world_up, forward).normalized();

        // if right length is zero (camera looking straight up/down), fix by choosing another up
        if (right.length() < 1e-6) {
            // use (0,1,0) as up instead
            world_up = Vec3(0,1,0);
            right = cross(world_up, forward).normalized();
        }

        Vec3 up = cross(forward, right);

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

    // Build a regex to find the camera block under CAMERA -> PERSP -> camera_id
    // This is a minimal, fragile JSON parser suitable for controlled test files
    // Pattern extracts numeric fields and nested vec objects with x,y,z
    const std::string idPattern = "\"" + camera_id + "\"\s*:\\s*\\{([\\s\\S]*?)\\}";
    std::regex id_regex(idPattern);
    std::smatch id_match;

    // First, locate the PERSP section
    std::regex persp_regex("\"CAMERA\"\s*:\\s*\\{[\\s\\S]*?\"PERSP\"\s*:\\s*\\{([\\s\\S]*?)\\}\\s*\\}");
    std::smatch persp_match;
    if (!std::regex_search(content, persp_match, persp_regex)) {
        return false;
    }
    const std::string persp_block = persp_match[1].str();

    if (!std::regex_search(persp_block, id_match, id_regex)) {
        return false; // camera id not found
    }
    const std::string cam_block = id_match[1].str();

    auto extract_vec3 = [](const std::string& src, const std::string& key, Vec3& out) -> bool {
        std::regex r("\"" + key + "\"\\s*:\\s*\\{\\s*\"x\"\\s*:\\s*([-+eE0-9\.]+)\\s*,\\s*\"y\"\\s*:\\s*([-+eE0-9\.]+)\\s*,\\s*\"z\"\\s*:\\s*([-+eE0-9\.]+)\\s*\\}");
        std::smatch m; if (!std::regex_search(src, m, r)) return false; out = Vec3(std::stod(m[1]), std::stod(m[2]), std::stod(m[3])); return true; };
    auto extract_number = [](const std::string& src, const std::string& key, double& out) -> bool {
        std::regex r("\"" + key + "\"\\s*:\\s*([-+eE0-9\.]+)");
        std::smatch m; if (!std::regex_search(src, m, r)) return false; out = std::stod(m[1]); return true; };
    auto extract_int = [](const std::string& src, const std::string& key, int& out) -> bool {
        std::regex r("\"" + key + "\"\\s*:\\s*([0-9]+)");
        std::smatch m; if (!std::regex_search(src, m, r)) return false; out = std::stoi(m[1]); return true; };
    auto extract_vec2i = [](const std::string& src, const std::string& key, int& outx, int& outy) -> bool {
        std::regex r("\"" + key + "\"\\s*:\\s*\\{\\s*\"x\"\\s*:\\s*([0-9]+)\\s*,\\s*\"y\"\\s*:\\s*([0-9]+)\\s*\\}");
        std::smatch m; if (!std::regex_search(src, m, r)) return false; outx = std::stoi(m[1]); outy = std::stoi(m[2]); return true; };

    Vec3 location, direction;
    double focal = 0.0, sensor_w = 0.0, sensor_h = 0.0; int resx = 0, resy = 0;
    bool ok = true;
    ok &= extract_vec3(cam_block, "location", location);
    ok &= extract_vec3(cam_block, "direction", direction);
    ok &= extract_number(cam_block, "focal_length", focal);
    ok &= extract_number(cam_block, "sensor_width", sensor_w);
    ok &= extract_number(cam_block, "sensor_height", sensor_h);
    ok &= extract_vec2i(cam_block, "film_resolution", resx, resy);
    if (!ok) return false;

    camera.setFocalLength(focal);
    camera.setSensorSize(sensor_w, sensor_h);
    camera.setImageResolution(resx, resy);
    camera.updateRotationMatrix(location, direction);
    return true;
}
