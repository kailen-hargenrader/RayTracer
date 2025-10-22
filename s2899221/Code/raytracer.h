#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include "mesh.h"
#include "utils.h"
#include "persp_camera.h"

/** Simple point light representation. */
struct PointLight {
    RayVec3 location;
    double radiant_intensity;
};

/** Scene container and core ray tracing utilities. */
class RayTracer {
public:
    RayTracer();

    /** Load meshes, cameras, and lights from a JSON file into this scene. */
    bool load_from_json(const std::string& json_filepath);

    /** Generate one primary ray per pixel for the camera with given id. */
    std::vector<Ray> get_one_ray_per_pixel(const std::string& camera_id) const;

    /** Intersect a single ray against all meshes; return the closest hit if any. */
    bool one_pass_intersection(const Ray& ray, Hit& out_hit) const;

    /** Access flat list of scene mesh pointers (for acceleration structures). */
    const std::vector<const Mesh*>& get_scene_mesh_ptrs() const { return m_scene_mesh_ptrs; }

private:
    // Store concrete meshes by value and maintain a flat list of pointers for polymorphic intersect
    std::vector<Cube> m_cubes;
    std::vector<Plane> m_planes;
    std::vector<Cylinder> m_cylinders;
    std::vector<Sphere> m_spheres;
    std::vector<const Mesh*> m_scene_mesh_ptrs;

    // Cameras addressable by id
    std::unordered_map<std::string, Camera> m_cameras;

    // Point lights addressable by id
    std::unordered_map<std::string, PointLight> m_point_lights;

    // Helper to rebuild pointer list after loading
    void rebuild_mesh_ptrs();
};

/**
 * Load a scene from JSON, generate one ray per pixel for the given camera id,
 * and time how long it takes to intersect all rays once. Returns seconds.
 */
double time_trace_all_one_pass(const std::string& json_filepath, const std::string& camera_id);

// -------- BVH acceleration API --------

class BoundingBox; // from mesh.h

/** Build a BVH using median splitting over object AABBs. */
const Mesh* get_BVH_tree(const std::vector<const Mesh*>& objects, std::vector<std::unique_ptr<BoundingBox>>& node_store);

/** Traverse BVH to find closest hit, returns true on any hit. */
bool accelerated_one_pass_intersection(const Mesh* bvh_root, const Ray& ray, Hit& out_hit);

/** Time ray tracing using the accelerated BVH traversal. Returns seconds. */
double time_trace_all_accelerated_one_pass(const std::string& json_filepath, const std::string& camera_id);

#endif // RAYTRACER_H
