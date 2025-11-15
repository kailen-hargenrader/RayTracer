## PPM Image Tests: Generate, Duplicate, Redshift, View

This folder contains small utilities to create a reference PPM image, duplicate it, apply a simple “redshift” effect, and visualize the results.

### 1) Generate a reference PPM
Use the Python helper to write a color map PPM:

```powershell
python tests\generate_colormap_ppm.py --output tests\ppm_image\reference_64x64.ppm
```

Notes:
- The generator writes an ASCII P3 PPM by default.
- The output file name `reference_64x64.ppm` is expected by the C++ drivers below.

---

### 2) Build the drivers with the Makefile
Build using the tests Makefile (outputs go to the repo root):

```powershell
cd tests
nmake ppm_image
```

---

### 3) Create the duplicate and redshifted images
Run the built executables (from repo root):

```powershell
.\duplicate_driver.exe
.\redshift_driver.exe
```

Outputs:
- `tests/ppm_image/duplicate_reference_64x64.ppm`
- `tests/ppm_image/redshift_reference_64x64.ppm`

Notes:
- The C++ `Image` writer outputs binary P6 PPM; the reference may be P3. This is expected—pixel values will match while file encodings differ.

---

### 4) Visualize the images
Use the PPM viewer script to quickly inspect images. If you have a viewer script (e.g., `ppm_viewer.py`), open it and point it at any of:
- `tests/ppm_image/reference_64x64.ppm`
- `tests/ppm_image/duplicate_reference_64x64.ppm`
- `tests/ppm_image/redshift_reference_64x64.ppm`

Example (PowerShell):

```powershell
python tests\raytracer\ppm_viewer.py
```

Or adjust the viewer script to load one of the files above and run it similarly.


