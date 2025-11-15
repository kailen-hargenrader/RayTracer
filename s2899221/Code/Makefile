CXX=cl
CXXFLAGS=/nologo /EHsc /std:c++17

# Build from s2899221/Code/ directory; write outputs to repo root (..\..)

all: ..\..\persp_camera.obj ..\..\mesh.obj ..\..\utils.obj ..\..\render_driver.obj ..\..\render_driver.exe

..\..\persp_camera.obj: persp_camera.cpp persp_camera.h utils.h
	$(CXX) $(CXXFLAGS) /c persp_camera.cpp /Fo:..\\..\\persp_camera.obj

..\..\mesh.obj: mesh.cpp mesh.h utils.h
	$(CXX) $(CXXFLAGS) /c mesh.cpp /Fo:..\\..\\mesh.obj

..\..\utils.obj: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) /c utils.cpp /Fo:..\\..\\utils.obj

..\..\render_driver.obj: render_driver.cpp persp_camera.h mesh.h utils.h raytracer.h
	$(CXX) $(CXXFLAGS) /c render_driver.cpp /Fo:..\\..\\render_driver.obj

..\..\render_driver.exe: render_driver.cpp
	$(CXX) $(CXXFLAGS) render_driver.cpp /Fe:..\\..\\render_driver.exe

clean:
	-del /Q ..\\..\\camera.obj ..\\..\\persp_camera.obj ..\\..\\mesh.obj ..\\..\\utils.obj ..\\..\\render_driver.obj 2>NUL
	-del /Q ..\\..\\render_driver.exe 2>NUL

