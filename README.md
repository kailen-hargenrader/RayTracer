# RayTracer

Building a ray tracer from scratch in C++.

This repository contains:
- A Blender export script to write scene data (cameras, meshes, lights) to JSON
- A perspective camera implementation to project world coordinates to pixels
- A PPM image utility for loading, modifying, and saving images

## Contents

- `s2899221/Blend/Export.py` — Blender script that exports scene data (CAMERA, MESH, LIGHT) to a JSON file.
- `s2899221/Code/persp_camera.h/.cpp` — Perspective camera and utilities to load camera parameters from JSON and compute world-to-pixel coordinates.
- `s2899221/Code/image.h/.cpp` — PPM image loader/writer with simple per-pixel manipulation.
- `tests/ppm_image` — Small utilities and drivers for PPM testing and visualization.

---

## 1) Export JSON from Blender using `Export.py`

Prerequisites:
- Blender installed

Steps:
1. Open your .blend scene in Blender.
2. Open the Scripting workspace (or Text Editor) in Blender.
3. Load `s2899221/Blend/Export.py` into Blender’s text editor.
4. Run the script inside Blender (Text Editor → Run Script).

What gets exported:
- Cameras (type PERSP): world `location`, normalized `direction`, `focal_length`, `sensor_width`, `sensor_height`, and `film_resolution` from render settings.
- Meshes: cubes (translation, rotation, scale), spheres (location, radius), planes (corner points).
- Lights: point lights with position and radiant intensity.

Where the file goes:
- By default, the script writes to: `~/RayTracer/s2899221/ASCII/cubes.json`.
- You can modify the target path at the bottom of `Export.py` (variable `filepath`).

---

## 2) Load camera from JSON and project world → pixel

Files:
- `s2899221/Code/persp_camera.h/.cpp`

Key APIs:
- `bool read_camera_from_json(const std::string& json_path, const std::string& camera_id, Camera& camera);`
- `bool Camera::worldToPixel(const Vec3& world, int& px, int& py);`

Usage example (MSVC, from repo root):

```bat
cd s2899221\Code
nmake
```

Minimal C++ snippet showing usage:

```cpp
#include <iostream>
#include "persp_camera.h"

int main() {
    Camera cam;
    // Use a camera id shown in your exported JSON (e.g., "1").
    if (!read_camera_from_json("s2899221/ASCII/cameras.json", "1", cam)) {
        std::cerr << "Failed to load camera" << std::endl;
        return 1;
    }

    Vec3 worldPoint{0.0, 0.0, 0.0};
    int px = 0, py = 0;
    if (cam.worldToPixel(worldPoint, px, py)) {
        std::cout << "Pixel: (" << px << ", " << py << ")\n";
    } else {
        std::cout << "Point is behind camera\n";
    }
}
```

Build (MSVC):

```bat
cl /nologo /EHsc /std:c++17 s2899221\Code\persp_camera.cpp s2899221\Code\camera.cpp /Fe:camera_demo.exe
```

Run the resulting executable as needed.

---

## 3) Load PPM image, modify pixels, write PPM

Files:
- `s2899221/Code/image.h/.cpp`

Capabilities:
- Read PPM P3 (ASCII) and P6 (binary) with max value in [1, 255].
- Store max value; clamp operations to 0..maxValue().
- Write PPM P6 using the current max value.

Basic usage:

```cpp
#include <iostream>
#include "image.h"

int main() {
    Image img("tests/ppm_image/reference_64x64.ppm");

    // Increase red channel (clamped)
    for (int y = 0; y < img.height(); ++y) {
        for (int x = 0; x < img.width(); ++x) {
            img.addToPixel(x, y, 32, 0, 0);
        }
    }

    img.write("tests/ppm_image/redshift_reference_64x64.ppm");
}
```

Build (MSVC):

```bat
cl /nologo /EHsc /std:c++17 s2899221\Code\image.cpp /Fe:image_demo.exe
```

Or use the provided Makefile (MSVC/nmake):

```bat
cd s2899221\Code
nmake
```

This produces `image.obj` and `persp_camera.obj` object files you can link into your own programs.

---

## 4) Test utilities (optional)

- `tests/ppm_image/generate_64x64_ppm.py` — writes a sample P3 PPM.
- `tests/ppm_image/ppm_viewer.py` — opens a PPM using Pillow. Default file is `reference_64x64.ppm` in the same folder.
- `tests/ppm_image/duplicate_driver.cpp` — loads a PPM with `Image` and writes a duplicate (P6).
- `tests/ppm_image/redshift_driver.cpp` — adds red to all pixels and writes a new PPM.

Example MSVC builds:

```bat
cl /nologo /EHsc /std:c++17 tests\ppm_image\duplicate_driver.cpp s2899221\Code\image.cpp /Fe:tests\ppm_image\duplicate_driver.exe
cl /nologo /EHsc /std:c++17 tests\ppm_image\redshift_driver.cpp s2899221\Code\image.cpp /Fe:tests\ppm_image\redshift_driver.exe
```

Run from repo root:

```bat
tests\ppm_image\duplicate_driver.exe
tests\ppm_image\redshift_driver.exe
```

Notes:
- Duplicate output is P6; original reference may be P3, so file bytes differ though pixel values correspond to the image's `maxValue()` range.

