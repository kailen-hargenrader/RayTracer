// driver.cpp
#include <iostream>
// Include implementation via relative path from this file's directory
#include "../../s2899221/Code/persp_camera.cpp"  // brings in Camera + read_camera_from_json

int main() {
    Camera cam;
    if (!read_camera_from_json("s2899221/ASCII/cameras.json", "1", cam)) {
        std::cerr << "Camera 1 not found\n";
        return 1;
    }

    Vec3 worldPoint(1.0, 2.0, 1.5); // replace with your point
    int px = 0, py = 0;
    if (cam.worldToPixel(worldPoint, px, py)) {
        std::cout << "Pixel: (" << px << ", " << py << ")\n";
    } else {
        std::cout << "Point is behind the camera\n";
    }
    return 0;
}