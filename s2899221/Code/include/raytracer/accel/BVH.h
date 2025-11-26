#pragma once

#include <memory>
#include <vector>
#include <cstdint>
#include "raytracer/accel/AABB.h"
#include "raytracer/utils/Hit.h"
#include "raytracer/utils/Ray.h"
#include "raytracer/geometry/Mesh.h"
#include "raytracer/core/Material.h"

namespace rt {

struct SceneObject {
    MeshPtr mesh;
    Material material;
    AABB bounds; // world-space
};

struct BVHNode {
    AABB bounds;
    int left {-1};
    int right {-1};
    int start {0};
    int count {0};
    bool isLeaf() const { return count > 0; }
};

class BVH {
public:
    // Build BVH over objects. Stores only indices internally.
    void build(const std::vector<SceneObject>& objects);

    // Intersect ray with BVH; returns true if any closer hit than tMax found and updates hit.
    bool intersect(const std::vector<SceneObject>& objects, const Ray& ray, float tMin, float tMax, Hit& hit, int* outObjectIndex = nullptr) const;

private:
    std::vector<BVHNode> m_nodes;
    std::vector<int> m_indices; // index into objects

    int buildRecursive(int start, int end, const std::vector<SceneObject>& objects, int depth);
};

} // namespace rt


