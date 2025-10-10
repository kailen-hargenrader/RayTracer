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

// Extracts the inner contents of an object block that follows a key in JSON-like text.
// Example: given src containing "\"FOO\": { ... }" and key "FOO",
// this returns the substring between the braces of that object (without the braces).
static bool extract_object_block(const std::string& src, const std::string& key, std::string& out_block) {
    const std::string quoted_key = "\"" + key + "\"";
    const size_t key_pos = src.find(quoted_key);
    if (key_pos == std::string::npos) return false;

    // Find the colon after the key
    size_t colon_pos = src.find(':', key_pos + quoted_key.size());
    if (colon_pos == std::string::npos) return false;

    // Find the opening brace of the object value
    size_t brace_start = src.find('{', colon_pos + 1);
    if (brace_start == std::string::npos) return false;

    // Walk forward and match braces, skipping over string contents
    int depth = 0;
    bool in_string = false;
    bool is_escaped = false;
    for (size_t i = brace_start; i < src.size(); ++i) {
        const char c = src[i];
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

        if (c == '"') {
            in_string = true;
            continue;
        }

        if (c == '{') {
            if (depth == 0) brace_start = i;
            ++depth;
        } else if (c == '}') {
            --depth;
            if (depth == 0) {
                const size_t brace_end = i;
                out_block = src.substr(brace_start + 1, brace_end - brace_start - 1);
                return true;
            }
        }
    }
    return false;
}

// Parses a Vec3 from a nested object block like { "x": <num>, "y": <num>, "z": <num> }
static bool parse_vec3_block(const std::string& src, const std::string& key, Vec3& out) {
    std::string block;
    if (!extract_object_block(src, key, block)) return false;

    std::smatch m;
    std::regex rx_x("\\\"x\\\"\\s*:\\s*([-+eE0-9\\.]+)");
    std::regex rx_y("\\\"y\\\"\\s*:\\s*([-+eE0-9\\.]+)");
    std::regex rx_z("\\\"z\\\"\\s*:\\s*([-+eE0-9\\.]+)");
    std::smatch mx, my, mz;
    if (!std::regex_search(block, mx, rx_x)) return false;
    if (!std::regex_search(block, my, rx_y)) return false;
    if (!std::regex_search(block, mz, rx_z)) return false;
    out = Vec3(std::stod(mx[1]), std::stod(my[1]), std::stod(mz[1]));
    return true;
}

// Parses two ints from a nested object block like { "x": <int>, "y": <int> }
static bool parse_vec2i_block(const std::string& src, const std::string& key, int& outx, int& outy) {
    std::string block;
    if (!extract_object_block(src, key, block)) return false;

    std::regex rx("\\\"x\\\"\\s*:\\s*([0-9]+)[^0-9]+\\\"y\\\"\\s*:\\s*([0-9]+)");
    std::smatch m;
    if (!std::regex_search(block, m, rx)) return false;
    outx = std::stoi(m[1]);
    outy = std::stoi(m[2]);
    return true;
}

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

    // Extract blocks step-by-step: CAMERA -> PERSP -> camera_id
    std::string camera_block;
    if (!extract_object_block(content, "CAMERA", camera_block)) {
        return false;
    }
    std::string persp_block;
    if (!extract_object_block(camera_block, "PERSP", persp_block)) {
        return false;
    }
    std::string cam_block;
    if (!extract_object_block(persp_block, camera_id, cam_block)) {
        return false;
    }

    auto extract_number = [](const std::string& src, const std::string& key, double& out) -> bool {
        std::regex r("\\\"" + key + "\\\"\\s*:\\s*([-+eE0-9\\.]+)");
        std::smatch m; if (!std::regex_search(src, m, r)) return false; out = std::stod(m[1]); return true; };

    Vec3 location, direction;
    double focal = 0.0, sensor_w = 0.0, sensor_h = 0.0; int resx = 0, resy = 0;
    bool ok = true;
    ok &= parse_vec3_block(cam_block, "location", location);
    ok &= parse_vec3_block(cam_block, "direction", direction);
    ok &= extract_number(cam_block, "focal_length", focal);
    ok &= extract_number(cam_block, "sensor_width", sensor_w);
    ok &= extract_number(cam_block, "sensor_height", sensor_h);
    ok &= parse_vec2i_block(cam_block, "film_resolution", resx, resy);
    if (!ok) return false;

    camera.setFocalLength(focal);
    camera.setSensorSize(sensor_w, sensor_h);
    camera.setImageResolution(resx, resy);
    camera.updateRotationMatrix(location, direction);
    return true;
}
