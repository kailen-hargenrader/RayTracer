cmake_minimum_required(VERSION 3.15)
project(RayTracerCore LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(raytracer STATIC
    Code/src/math/Mat4.cpp
    Code/src/geometry/Sphere.cpp
    Code/src/geometry/Cube.cpp
    Code/src/geometry/Cylinder.cpp
    Code/src/geometry/Plane.cpp
    Code/src/utils/Image.cpp
    Code/src/utils/Texture.cpp
    Code/src/utils/Json.cpp
    Code/src/accel/AABB.cpp
    Code/src/accel/BVH.cpp
    Code/src/core/Scene.cpp
    Code/src/core/Renderer.cpp
)

target_include_directories(raytracer PUBLIC Code/include)

add_executable(rt_render
    Code/apps/rt_render.cpp
)
target_link_libraries(rt_render PRIVATE raytracer)


