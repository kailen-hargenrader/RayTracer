## Intersections Test: Generate, Load, and Visualize

This folder contains tools to:
- Generate random camera rays from the first perspective camera in a scene JSON and compute their intersections with scene geometry (`test_intersect_fn.cpp`).
- Load the resulting ray segments into Blender as cylinders (`load_intersections.py`).
- Visualize intersections directly inside Blender from a scene JSON (`visualize_intersections.py`).

### Prerequisites
- A scene JSON in the project (e.g., files under `s2899221/ASCII/`).
- C++17 compiler (for building `test_intersect_fn`).
- Blender installed (accessible from your Shell or via full path).

---

### 1) Build and run `test_intersect_fn` (generate ray intersections)
This C++ program:
- Reads the first perspective camera from the input JSON.
- Generates random pixel rays using that camera’s film resolution and intrinsics.
- Intersects those rays with all meshes in the scene (Cubes, Planes, Cylinders, Spheres).
- Writes out line segments representing the primary ray to the first hit (and a reflected segment if applicable).

Build with the provided Makefile (outputs go to repo root):

```powershell
cd tests
nmake intersect
```

Usage:

```powershell
.\test_intersect_fn.exe <path-to-scene.json> <num-rays> [out.txt] [max_length]
```

- `<path-to-scene.json>`: Input scene JSON.
- `<num-rays>`: How many random rays to sample from the first camera.
- `[out.txt]` (optional): Output file for segments; defaults to `intersections.txt` in the working directory.
- `[max_length]` (optional): Maximum length to cap rays; defaults to `8.0`.

Output format (each line is one cylinder segment):
```
ox oy oz dx dy dz L
```
- `origin = (ox, oy, oz)`
- `direction = (dx, dy, dz)` (unit or a delta depending on context)
- `length = L`

Notes:
- On a hit within `max_length`, the tool writes:
  - One segment from the ray origin to the hit point
  - One reflected segment starting at the hit point (remaining length)
- If no hit within `max_length`, it writes a single segment for the primary ray capped at `max_length`.

Example run:

```powershell
.\test_intersect_fn.exe s2899221\ASCII\lab2.json 300 intersections.txt 8.0
```

---

### 2) Load intersections into Blender (`load_intersections.py`)
This script reads the segment file produced by `test_intersect_fn` and creates cylinders along each segment in Blender.

Before running:
- Open `tests/intersect/load_intersections.py`
- Set `FILE_PATH` to the segment file you generated (e.g., `C:\\Users\\Kailen\\RayTracer\\intersections.txt`)
- Optionally adjust:
  - `CYL_RADIUS` (cylinder radius)
  - `CREATE_COLLECTION` and `COLLECTION_NAME`

Run in Blender (PowerShell example):

```powershell
blender --factory-startup --python tests\intersect\load_intersections.py
```

Result:
- Cylinders are created for each segment in the specified collection. You should see primary and reflected segments where applicable.

---

### 3) Visualize intersections directly in Blender (`visualize_intersections.py`)
This script computes intersections inside Blender starting from a scene JSON. It:
- Reads the first perspective camera from the JSON
- Generates random rays and intersects them with supported meshes (Cubes, Planes)
- Draws cylinders for primary rays and their short reflection rays

Before running:
- Open `tests/intersect/visualize_intersections.py`
- Set `JSON_PATH` to your scene JSON (e.g., `C:\\Users\\Kailen\\RayTracer\\s2899221\\ASCII\\lab2.json`)
- Optionally adjust:
  - `NUM_RAYS`, `MAX_RAY_LENGTH`, `REFLECT_LENGTH`
  - Cylinder radii and collection names

Run in Blender:

```powershell
blender --factory-startup --python tests\intersect\visualize_intersections.py
```

Verification:
- The ray/segment visualization from `load_intersections.py` (using the C++ output) and from `visualize_intersections.py` (computed in Blender) should look qualitatively similar for the same scene.
- Differences can arise from random sampling, supported mesh types, numeric tolerances, or parameter choices (`MAX_RAY_LENGTH`, `REFLECT_LENGTH`, radii).

---

### Tips
- Keep paths consistent; on Windows PowerShell, absolute paths with `\\` work best in Blender scripts.
- If your segment file is not in the repo root, update `FILE_PATH` accordingly.
- For reproducibility, keep `num-rays` and script parameters consistent across runs.


