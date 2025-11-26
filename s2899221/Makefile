!IFNDEF CC
CC=cl
!ENDIF

CXXFLAGS=/nologo /std:c++17 /EHsc /O2 /MD /I Code\include
LDFLAGS=/nologo

OUTDIR=build\bin
OBJDIR=build\obj

OBJECTS= \
	$(OBJDIR)\Mat4.obj \
	$(OBJDIR)\Sphere.obj \
	$(OBJDIR)\Cube.obj \
	$(OBJDIR)\Cylinder.obj \
	$(OBJDIR)\Plane.obj \
	$(OBJDIR)\Image.obj \
	$(OBJDIR)\Texture.obj \
	$(OBJDIR)\Json.obj \
	$(OBJDIR)\AABB.obj \
	$(OBJDIR)\BVH.obj \
	$(OBJDIR)\Scene.obj \
	$(OBJDIR)\Renderer.obj \
	$(OBJDIR)\rt_render.obj

all: dirs $(OUTDIR)\rt_render.exe

dirs:
	@if not exist "$(OUTDIR)" mkdir "$(OUTDIR)"
	@if not exist "$(OBJDIR)" mkdir "$(OBJDIR)"

$(OUTDIR)\rt_render.exe: $(OBJECTS)
	link $(LDFLAGS) /OUT:$(OUTDIR)\rt_render.exe $(OBJECTS)

$(OBJDIR)\Mat4.obj: Code\src\math\Mat4.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Mat4.obj Code\src\math\Mat4.cpp

$(OBJDIR)\Sphere.obj: Code\src\geometry\Sphere.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Sphere.obj Code\src\geometry\Sphere.cpp

$(OBJDIR)\Cube.obj: Code\src\geometry\Cube.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Cube.obj Code\src\geometry\Cube.cpp

$(OBJDIR)\Cylinder.obj: Code\src\geometry\Cylinder.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Cylinder.obj Code\src\geometry\Cylinder.cpp

$(OBJDIR)\Plane.obj: Code\src\geometry\Plane.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Plane.obj Code\src\geometry\Plane.cpp

$(OBJDIR)\Image.obj: Code\src\utils\Image.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Image.obj Code\src\utils\Image.cpp

$(OBJDIR)\Texture.obj: Code\src\utils\Texture.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Texture.obj Code\src\utils\Texture.cpp

$(OBJDIR)\Json.obj: Code\src\utils\Json.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Json.obj Code\src\utils\Json.cpp

$(OBJDIR)\AABB.obj: Code\src\accel\AABB.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\AABB.obj Code\src\accel\AABB.cpp

$(OBJDIR)\BVH.obj: Code\src\accel\BVH.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\BVH.obj Code\src\accel\BVH.cpp

$(OBJDIR)\Scene.obj: Code\src\core\Scene.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Scene.obj Code\src\core\Scene.cpp

$(OBJDIR)\Renderer.obj: Code\src\core\Renderer.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\Renderer.obj Code\src\core\Renderer.cpp

$(OBJDIR)\rt_render.obj: Code\apps\rt_render.cpp
	$(CC) /c $(CXXFLAGS) /Fo$(OBJDIR)\rt_render.obj Code\apps\rt_render.cpp

clean:
	@if exist "$(OBJDIR)" rmdir /S /Q "$(OBJDIR)"
	@if exist "$(OUTDIR)\rt_render.exe" del /Q "$(OUTDIR)\rt_render.exe"
	@echo Clean complete.


