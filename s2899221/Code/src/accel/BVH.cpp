#include "raytracer/accel/BVH.h"
#include <algorithm>

namespace rt {

static AABB unionOf(const AABB& a, const AABB& b) {
    AABB out = a;
    out.expand(b);
    return out;
}

void BVH::build(const std::vector<SceneObject>& objects) {
    m_indices.resize(objects.size());
    for (size_t i = 0; i < objects.size(); ++i) m_indices[i] = static_cast<int>(i);
    m_nodes.clear();
    if (!objects.empty()) {
        buildRecursive(0, static_cast<int>(objects.size()), objects, 0);
    }
}

int BVH::buildRecursive(int start, int end, const std::vector<SceneObject>& objects, int depth) {
    BVHNode node;
    AABB bounds;
    for (int i = start; i < end; ++i) {
        bounds.expand(objects[m_indices[i]].bounds);
    }
    node.bounds = bounds;
    int count = end - start;
    int nodeIndex = static_cast<int>(m_nodes.size());
    m_nodes.push_back(node);

    if (count <= 4 || depth > 32) {
        // Make leaf
        m_nodes[nodeIndex].start = start;
        m_nodes[nodeIndex].count = count;
        return nodeIndex;
    }

    // Choose split axis by extent
    Vec3 diag = bounds.max - bounds.min;
    int axis = 0;
    if (diag.y > diag.x && diag.y >= diag.z) axis = 1;
    else if (diag.z > diag.x) axis = 2;

    float mid = 0.5f * (bounds.min[axis] + bounds.max[axis]);
    int pivot = std::partition(m_indices.begin() + start, m_indices.begin() + end,
        [&](int idx) {
            const AABB& b = objects[idx].bounds;
            float center = 0.5f * (b.min[axis] + b.max[axis]);
            return center < mid;
        }) - m_indices.begin();
    if (pivot == start || pivot == end) {
        pivot = start + count / 2; // fallback median
    }

    int left = buildRecursive(start, pivot, objects, depth + 1);
    int right = buildRecursive(pivot, end, objects, depth + 1);

    m_nodes[nodeIndex].left = left;
    m_nodes[nodeIndex].right = right;
    return nodeIndex;
}

bool BVH::intersect(const std::vector<SceneObject>& objects, const Ray& ray, float tMin, float tMax, Hit& hit, int* outObjectIndex) const {
    if (m_nodes.empty()) return false;
    bool any = false;
    int bestIdx = -1;
    struct StackItem { int node; float tmin; float tmax; };
    StackItem stack[64];
    int sp = 0;
    stack[sp++] = {0, tMin, tMax};
    while (sp) {
        StackItem item = stack[--sp];
        const BVHNode& n = m_nodes[item.node];
        if (!n.bounds.intersect(ray, item.tmin, item.tmax)) continue;
        if (n.isLeaf()) {
            for (int i = 0; i < n.count; ++i) {
                int idx = m_indices[n.start + i];
                Hit tmp = hit;
                if (objects[idx].mesh->intersect(ray, item.tmin, hit.didHit ? hit.t : item.tmax, tmp)) {
                    // Carry material
                    tmp.albedo = objects[idx].material.baseColor;
                    hit = tmp;
                    bestIdx = idx;
                    any = true;
                }
            }
        } else {
            // Traverse children; push both
            stack[sp++] = { n.left, item.tmin, item.tmax };
            stack[sp++] = { n.right, item.tmin, item.tmax };
        }
    }
    if (any && outObjectIndex) *outObjectIndex = bestIdx;
    return any;
}

} // namespace rt


