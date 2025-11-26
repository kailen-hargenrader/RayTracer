#pragma once

#include <memory>
#include <vector>
#include <string>
#include "raytracer/core/Camera.h"
#include "raytracer/accel/BVH.h"
#include "raytracer/utils/Json.h"

namespace rt {

struct PointLight {
    Vec3 position;
    float intensity; // radiant intensity proxy
};

class Scene {
public:
    Camera camera;
    std::vector<PointLight> lights;
    std::vector<SceneObject> objects;
    BVH bvh; // optional acceleration

    static std::shared_ptr<Scene> loadFromJsonFile(const std::string& path);
};

} // namespace rt


