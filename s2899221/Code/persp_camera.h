#ifndef CAMERA_H
#define CAMERA_H

#include <string>

/** Simple 3D vector used by the camera math. */
struct Vec3 {
    double x, y, z;

    Vec3();
    Vec3(double _x, double _y, double _z);

    Vec3 operator-(const Vec3& v) const;
    Vec3 operator+(const Vec3& v) const;
    Vec3 operator*(double s) const;
    double length() const;
    Vec3 normalized() const;
};

/**
 * Perspective camera with simple pinhole projection.
 *
 * Fields include pose (position and rotation), lens parameters (focal length
 * and sensor size), and the output image resolution in pixels. Provides
 * world-to-pixel projection and a helper to read parameters from a minimal
 * JSON description used by the tests.
 */
class Camera {
private:
    Vec3 position;
    double rotation[3][3];  // columns: right, up, forward
    double focal_length_mm;
    double sensor_width_mm;
    double sensor_height_mm;
    int image_width_px;
    int image_height_px;

    Vec3 cross(const Vec3& a, const Vec3& b);

public:
    Camera();

    void setPosition(const Vec3& pos);
    void setFocalLength(double focal_mm);
    void setSensorSize(double width_mm, double height_mm);
    void setImageResolution(int width_px, int height_px);

    /** Update orientation from camera location and forward direction. */
    void updateRotationMatrix(const Vec3& cam_pos, const Vec3& cam_dir);

    /** Project a world point to pixel coordinates. Returns false if behind camera. */
    bool worldToPixel(const Vec3& point_world, int& pixel_x, int& pixel_y);
};

// Reads camera parameters from a JSON file and updates the provided Camera instance.
// Returns true if a camera with the given id was found and applied; false otherwise.
/** Parse and load a camera by id from the JSON file into the provided camera. */
bool read_camera_from_json(const std::string& json_filepath, const std::string& camera_id, Camera& camera);

#endif // CAMERA_H
