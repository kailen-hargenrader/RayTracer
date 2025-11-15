CXX=cl
CXXFLAGS=/nologo /EHsc /std:c++17

# Build from tests/ directory; write all .obj and .exe outputs to repo root (..)

# ----- Grouped by folder (explicit files) -----
# ppm_image
PPM_IMAGE_SOURCES=ppm_image\\duplicate_driver.cpp ppm_image\\redshift_driver.cpp
PPM_IMAGE_OBJECTS=..\\duplicate_driver.obj ..\\redshift_driver.obj
PPM_IMAGE_EXES=..\\duplicate_driver.exe ..\\redshift_driver.exe

# intersect
INTERSECT_SOURCES=intersect\\test_intersect_fn.cpp
INTERSECT_OBJECTS=..\\test_intersect_fn.obj
INTERSECT_EXES=..\\test_intersect_fn.exe

# raytracer
RAYTRACER_SOURCES=raytracer\\time_driver.cpp
RAYTRACER_OBJECTS=..\\time_driver.obj
RAYTRACER_EXES=

# core (from s2899221/Code)
CORE_OBJECTS=..\\utils.obj ..\\mesh.obj

.PHONY: all clean ppm_image intersect raytracer

all: ppm_image intersect raytracer

# ----- ppm_image (explicit) -----
ppm_image: $(PPM_IMAGE_EXES)

..\\duplicate_driver.obj: ppm_image\\duplicate_driver.cpp
	$(CXX) $(CXXFLAGS) /c ppm_image\\duplicate_driver.cpp /Fo:..\\duplicate_driver.obj

..\\redshift_driver.obj: ppm_image\\redshift_driver.cpp
	$(CXX) $(CXXFLAGS) /c ppm_image\\redshift_driver.cpp /Fo:..\\redshift_driver.obj

..\\duplicate_driver.exe: ..\\duplicate_driver.obj ..\\utils.obj
	$(CXX) $(CXXFLAGS) ..\\duplicate_driver.obj ..\\utils.obj /Fe:..\\duplicate_driver.exe

..\\redshift_driver.exe: ..\\redshift_driver.obj ..\\utils.obj
	$(CXX) $(CXXFLAGS) ..\\redshift_driver.obj ..\\utils.obj /Fe:..\\redshift_driver.exe

# ----- intersect (explicit) -----
intersect: $(INTERSECT_EXES)

..\\test_intersect_fn.obj: intersect\\test_intersect_fn.cpp
	$(CXX) $(CXXFLAGS) /c intersect\\test_intersect_fn.cpp /Fo:..\\test_intersect_fn.obj

..\\test_intersect_fn.exe: ..\\test_intersect_fn.obj ..\\mesh.obj ..\\utils.obj
	$(CXX) $(CXXFLAGS) ..\\test_intersect_fn.obj ..\\mesh.obj ..\\utils.obj /Fe:..\\test_intersect_fn.exe

# ----- raytracer (explicit) -----
raytracer: $(RAYTRACER_OBJECTS)

..\\time_driver.obj: raytracer\\time_driver.cpp
	$(CXX) $(CXXFLAGS) /c raytracer\\time_driver.cpp /Fo:..\\time_driver.obj

# ----- core objects (explicit) -----
..\\utils.obj: ..\\s2899221\\Code\\utils.cpp ..\\s2899221\\Code\\utils.h
	$(CXX) $(CXXFLAGS) /c ..\\s2899221\\Code\\utils.cpp /Fo:..\\utils.obj

..\\mesh.obj: ..\\s2899221\\Code\\mesh.cpp ..\\s2899221\\Code\\mesh.h ..\\s2899221\\Code\\utils.h
	$(CXX) $(CXXFLAGS) /c ..\\s2899221\\Code\\mesh.cpp /Fo:..\\mesh.obj

# ----- clean (explicit) -----
clean:
	-del /Q ..\\duplicate_driver.obj ..\\redshift_driver.obj 2>NUL
	-del /Q ..\\test_intersect_fn.obj 2>NUL
	-del /Q ..\\time_driver.obj 2>NUL
	-del /Q ..\\utils.obj ..\\mesh.obj 2>NUL
	-del /Q ..\\duplicate_driver.exe ..\\redshift_driver.exe ..\\test_intersect_fn.exe 2>NUL


