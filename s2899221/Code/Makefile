CXX=cl
CXXFLAGS=/nologo /EHsc /std:c++17

all: persp_camera.obj mesh.obj utils.obj

persp_camera.obj: persp_camera.cpp persp_camera.h utils.h
	$(CXX) $(CXXFLAGS) /c persp_camera.cpp /Fo:persp_camera.obj

mesh.obj: mesh.cpp mesh.h utils.h
	$(CXX) $(CXXFLAGS) /c mesh.cpp /Fo:mesh.obj

utils.obj: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) /c utils.cpp /Fo:utils.obj

clean:
	-del /Q camera.obj persp_camera.obj mesh.obj utils.obj 2>NUL

