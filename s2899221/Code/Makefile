CXX=cl
CXXFLAGS=/nologo /EHsc /std:c++17

all: image.obj persp_camera.obj mesh.obj

image.obj: image.cpp image.h
	$(CXX) $(CXXFLAGS) /c image.cpp /Fo:image.obj

persp_camera.obj: persp_camera.cpp persp_camera.h
	$(CXX) $(CXXFLAGS) /c persp_camera.cpp /Fo:persp_camera.obj

mesh.obj: mesh.cpp mesh.h
	$(CXX) $(CXXFLAGS) /c mesh.cpp /Fo:mesh.obj

clean:
	-del /Q image.obj camera.obj persp_camera.obj mesh.obj 2>NUL

