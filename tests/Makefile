CXX=cl
CXXFLAGS=/nologo /EHsc /std:c++17

# Explicitly listed sources by subdirectory (no recursion, nmake-friendly paths)
PPM_IMAGE_SOURCES=ppm_image\\duplicate_driver.cpp ppm_image\\redshift_driver.cpp
PERSP_CAMERA_SOURCES=persp_camera\\camera_to_rays.cpp persp_camera\\driver.cpp
INTERSECT_SOURCES=intersect\\test_intersect_fn.cpp
MESH_SOURCES=mesh\\json_to_console.cpp
RAYTRACER_SOURCES=raytracer\\time_driver.cpp raytracer\\render_driver.cpp

PPM_IMAGE_OBJECTS=ppm_image\\duplicate_driver.obj ppm_image\\redshift_driver.obj
PERSP_CAMERA_OBJECTS=persp_camera\\camera_to_rays.obj persp_camera\\driver.obj
INTERSECT_OBJECTS=intersect\\test_intersect_fn.obj
MESH_OBJECTS=mesh\\json_to_console.obj
RAYTRACER_OBJECTS=raytracer\\time_driver.obj raytracer\\render_driver.obj

ALL_OBJECTS=$(PPM_IMAGE_OBJECTS) $(PERSP_CAMERA_OBJECTS) $(INTERSECT_OBJECTS) $(MESH_OBJECTS) $(RAYTRACER_OBJECTS)

.PHONY: all clean ppm_image persp_camera intersect mesh raytracer

all: ppm_image persp_camera intersect mesh raytracer

# Per-subdirectory builds
ppm_image: $(PPM_IMAGE_OBJECTS)
persp_camera: $(PERSP_CAMERA_OBJECTS)
intersect: $(INTERSECT_OBJECTS)
mesh: $(MESH_OBJECTS)
raytracer: $(RAYTRACER_OBJECTS)
	$(CXX) $(CXXFLAGS) raytracer\\render_driver.obj /Fe:raytracer\\render_driver.exe

# Explicit compile rules for each file (avoids GNU make pattern rules)
ppm_image\\duplicate_driver.obj: ppm_image\\duplicate_driver.cpp
	$(CXX) $(CXXFLAGS) /c ppm_image\\duplicate_driver.cpp /Fo:ppm_image\\duplicate_driver.obj

ppm_image\\redshift_driver.obj: ppm_image\\redshift_driver.cpp
	$(CXX) $(CXXFLAGS) /c ppm_image\\redshift_driver.cpp /Fo:ppm_image\\redshift_driver.obj

persp_camera\\camera_to_rays.obj: persp_camera\\camera_to_rays.cpp
	$(CXX) $(CXXFLAGS) /c persp_camera\\camera_to_rays.cpp /Fo:persp_camera\\camera_to_rays.obj

persp_camera\\driver.obj: persp_camera\\driver.cpp
	$(CXX) $(CXXFLAGS) /c persp_camera\\driver.cpp /Fo:persp_camera\\driver.obj

intersect\\test_intersect_fn.obj: intersect\\test_intersect_fn.cpp
	$(CXX) $(CXXFLAGS) /c intersect\\test_intersect_fn.cpp /Fo:intersect\\test_intersect_fn.obj

mesh\\json_to_console.obj: mesh\\json_to_console.cpp
	$(CXX) $(CXXFLAGS) /c mesh\\json_to_console.cpp /Fo:mesh\\json_to_console.obj

raytracer\\time_driver.obj: raytracer\\time_driver.cpp
	$(CXX) $(CXXFLAGS) /c raytracer\\time_driver.cpp /Fo:raytracer\\time_driver.obj

raytracer\\render_driver.obj: raytracer\\render_driver.cpp
	$(CXX) $(CXXFLAGS) /c raytracer\\render_driver.cpp /Fo:raytracer\\render_driver.obj

raytracer\\render_driver.exe: raytracer\\render_driver.obj
	$(CXX) $(CXXFLAGS) raytracer\\render_driver.obj /Fe:raytracer\\render_driver.exe

clean:
	-del /Q $(ALL_OBJECTS) raytracer\\render_driver.exe 2>NUL


