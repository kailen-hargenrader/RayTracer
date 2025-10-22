#include "raytracer.h"

#include <fstream>
#include <sstream>
#include <regex>
#include <limits>
#include <chrono>
#include <algorithm>
#include <stack>

using util_json::extract_object_block;
using util_json::parse_number;
using util_json::parse_vec3;
using util_json::parse_vec2i;

RayTracer::RayTracer() = default;

void RayTracer::rebuild_mesh_ptrs() {
	m_scene_mesh_ptrs.clear();
	m_scene_mesh_ptrs.reserve(m_cubes.size() + m_planes.size() + m_cylinders.size() + m_spheres.size());
	for (const auto& c : m_cubes) m_scene_mesh_ptrs.push_back(&c);
	for (const auto& p : m_planes) m_scene_mesh_ptrs.push_back(&p);
	for (const auto& cy : m_cylinders) m_scene_mesh_ptrs.push_back(&cy);
	for (const auto& s : m_spheres) m_scene_mesh_ptrs.push_back(&s);
}

static bool read_all_cameras(const std::string& content, std::unordered_map<std::string, Camera>& out) {
	std::string camera_block; if (!extract_object_block(content, "CAMERA", camera_block)) return false;
	std::string persp_block; if (!extract_object_block(camera_block, "PERSP", persp_block)) return false;

	// Iterate over all id objects in PERSP block
	std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
	auto begin = std::sregex_iterator(persp_block.begin(), persp_block.end(), id_rx);
	auto end = std::sregex_iterator();
	for (auto it = begin; it != end; ++it) {
		const size_t key_pos = static_cast<size_t>(it->position());
		std::string tail = persp_block.substr(key_pos);
		std::string cam_block; if (!extract_object_block(tail, it->str(1), cam_block)) continue;

		Camera cam;
		// Reuse camera parsing fields
		double x=0,y=0,z=0;
		if (!parse_vec3(cam_block, "location", x, y, z)) continue;
		Vec3 location(x,y,z);
		if (!parse_vec3(cam_block, "direction", x, y, z)) continue;
		Vec3 direction(x,y,z);
		double focal=35.0, sw=36.0, sh=24.0; int rx=1920, ry=1080;
		parse_number(cam_block, "focal_length", focal);
		parse_number(cam_block, "sensor_width", sw);
		parse_number(cam_block, "sensor_height", sh);
		parse_vec2i(cam_block, "film_resolution", rx, ry);
		cam.setFocalLength(focal);
		cam.setSensorSize(sw, sh);
		cam.setImageResolution(rx, ry);
		cam.updateRotationMatrix(location, direction);
		out[it->str(1)] = cam;
	}
	return !out.empty();
}

static void read_all_meshes(const std::string& content, std::vector<Cube>& cubes, std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres) {
	std::string mesh_block; if (!extract_object_block(content, "MESH", mesh_block)) return;
	std::string cube_block; if (extract_object_block(mesh_block, "Cube", cube_block)) {
		auto v = Cube::read_from_json(cube_block); cubes.insert(cubes.end(), v.begin(), v.end());
	}
	std::string plane_block; if (extract_object_block(mesh_block, "Plane", plane_block)) {
		auto v = Plane::read_from_json(plane_block); planes.insert(planes.end(), v.begin(), v.end());
	}
	std::string cyl_block; if (extract_object_block(mesh_block, "Cylinder", cyl_block)) {
		auto v = Cylinder::read_from_json(cyl_block); cylinders.insert(cylinders.end(), v.begin(), v.end());
	}
	std::string sph_block; if (extract_object_block(mesh_block, "Sphere", sph_block)) {
		auto v = Sphere::read_from_json(sph_block); spheres.insert(spheres.end(), v.begin(), v.end());
	}
}

static void read_point_lights(const std::string& content, std::unordered_map<std::string, PointLight>& out) {
	std::string light_block; if (!extract_object_block(content, "LIGHT", light_block)) return;
	std::string point_block; if (!extract_object_block(light_block, "POINT", point_block)) return;
	std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
	auto begin = std::sregex_iterator(point_block.begin(), point_block.end(), id_rx);
	auto end = std::sregex_iterator();
	for (auto it = begin; it != end; ++it) {
		const size_t key_pos = static_cast<size_t>(it->position());
		std::string tail = point_block.substr(key_pos);
		std::string pl_block; if (!extract_object_block(tail, it->str(1), pl_block)) continue;
		double x=0,y=0,z=0; if (!parse_vec3(pl_block, "location", x, y, z)) continue;
		double I = 1.0; parse_number(pl_block, "radiant_intensity", I);
		PointLight L; L.location = {x,y,z}; L.radiant_intensity = I;
		out[it->str(1)] = L;
	}
}

bool RayTracer::load_from_json(const std::string& json_filepath) {
	std::ifstream in(json_filepath);
	if (!in) return false;
	std::stringstream buffer; buffer << in.rdbuf();
	const std::string content = buffer.str();

	m_cubes.clear(); m_planes.clear(); m_cylinders.clear(); m_spheres.clear();
	m_cameras.clear(); m_point_lights.clear(); m_scene_mesh_ptrs.clear();

	read_all_meshes(content, m_cubes, m_planes, m_cylinders, m_spheres);
	read_all_cameras(content, m_cameras);
	read_point_lights(content, m_point_lights);
	rebuild_mesh_ptrs();
	return true;
}

std::vector<Ray> RayTracer::get_one_ray_per_pixel(const std::string& camera_id) const {
	std::vector<Ray> rays;
	auto it = m_cameras.find(camera_id);
	if (it == m_cameras.end()) return rays;
	const Camera& cam = it->second;
	int width=0, height=0; cam.getImageResolution(width, height);
	rays.reserve(static_cast<size_t>(width) * static_cast<size_t>(height));

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			Vec3 ro, rd;
			cam.pixelToRay(static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5, ro, rd);
			RayVec3 o{ ro.x, ro.y, ro.z };
			RayVec3 d{ rd.x, rd.y, rd.z };
			rays.emplace_back(o, d);
		}
	}
	return rays;
}

bool RayTracer::one_pass_intersection(const Ray& ray, Hit& out_hit) const {
	bool any_hit = false;
	double best_dist = std::numeric_limits<double>::infinity();
	Hit best_hit;
	for (const Mesh* m : m_scene_mesh_ptrs) {
		Hit h;
		if (m->intersect(ray, h) && h.hasDistanceAlongRay()) {
			double dist = h.getDistanceAlongRay();
			if (dist < best_dist) { best_dist = dist; best_hit = h; any_hit = true; }
		}
	}
	if (any_hit) { out_hit = best_hit; }
	return any_hit;
}

// ---------------- BVH construction and traversal ----------------

static inline Float3 aabb_union_min(const Float3& a, const Float3& b) { return Float3{ std::min(a.x,b.x), std::min(a.y,b.y), std::min(a.z,b.z) }; }
static inline Float3 aabb_union_max(const Float3& a, const Float3& b) { return Float3{ std::max(a.x,b.x), std::max(a.y,b.y), std::max(a.z,b.z) }; }

static inline double centroid_axis(const Float3& mn, const Float3& mx, int axis) {
    const double c[3] = { 0.5*(mn.x+mx.x), 0.5*(mn.y+mx.y), 0.5*(mn.z+mx.z) };
    return c[axis];
}

const Mesh* get_BVH_tree(const std::vector<const Mesh*>& objects, std::vector<std::unique_ptr<BoundingBox>>& node_store) {
    if (objects.empty()) return nullptr;
    if (objects.size() == 1) return objects[0];

    // Compute bounds of all objects and their centroids
    std::vector<std::pair<const Mesh*, std::pair<Float3,Float3>>> items;
    items.reserve(objects.size());
    Float3 scene_min{ +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity() };
    Float3 scene_max{ -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
    for (const Mesh* m : objects) {
        Float3 mn, mx; m->compute_aabb(mn, mx);
        items.emplace_back(m, std::make_pair(mn, mx));
        scene_min = aabb_union_min(scene_min, mn);
        scene_max = aabb_union_max(scene_max, mx);
    }

    // Choose split axis by largest extent
    const double extents[3] = { scene_max.x - scene_min.x, scene_max.y - scene_min.y, scene_max.z - scene_min.z };
    int axis = 0; if (extents[1] > extents[axis]) axis = 1; if (extents[2] > extents[axis]) axis = 2;

    // Sort by centroid along axis and split by median
    std::sort(items.begin(), items.end(), [&](const auto& A, const auto& B){
        return centroid_axis(A.second.first, A.second.second, axis) < centroid_axis(B.second.first, B.second.second, axis);
    });
    const size_t mid = items.size()/2;
    std::vector<const Mesh*> left_objs; left_objs.reserve(mid);
    std::vector<const Mesh*> right_objs; right_objs.reserve(items.size()-mid);
    for (size_t i = 0; i < items.size(); ++i) {
        if (i < mid) left_objs.push_back(items[i].first); else right_objs.push_back(items[i].first);
    }

    // Recursively build children
    const Mesh* left = get_BVH_tree(left_objs, node_store);
    const Mesh* right = get_BVH_tree(right_objs, node_store);

    // Create internal node bounding both children
    Float3 lmin{0,0,0}, lmax{0,0,0}, rmin{0,0,0}, rmax{0,0,0};
    left->compute_aabb(lmin, lmax);
    right->compute_aabb(rmin, rmax);
    Float3 nmin = aabb_union_min(lmin, rmin);
    Float3 nmax = aabb_union_max(lmax, rmax);
    node_store.emplace_back(new BoundingBox(nmin, nmax));
    BoundingBox* node = node_store.back().get();
    node->setChildren(left, right);
    return node;
}

bool accelerated_one_pass_intersection(const Mesh* bvh_root, const Ray& ray, Hit& out_hit) {
    if (!bvh_root) return false;
    // Iterative stack for traversal; rely on BoundingBox::intersect to test AABB quickly
    struct StackItem { const Mesh* node; double tnear; };
    std::stack<StackItem> st;
    st.push({ bvh_root, 0.0 });
    bool any_hit = false;
    double best_dist = std::numeric_limits<double>::infinity();
    Hit best_hit;

    while (!st.empty()) {
        const Mesh* node = st.top().node; st.pop();
        auto bbox_ptr = node->asBoundingBox();
        if (bbox_ptr) {
            Hit aabb_hit;
            if (!bbox_ptr->intersect(ray, aabb_hit)) continue; // prune
            if (aabb_hit.hasDistanceAlongRay() && aabb_hit.getDistanceAlongRay() > best_dist) continue; // farther than best hit
            const Mesh* left = bbox_ptr->getLeft();
            const Mesh* right = bbox_ptr->getRight();
            if (left && right) {
                // Order children by their AABB entry distance to visit nearer first
                Float3 lmin{0,0,0}, lmax{0,0,0}, rmin{0,0,0}, rmax{0,0,0};
                left->compute_aabb(lmin, lmax); right->compute_aabb(rmin, rmax);
                // Approximate tnear by intersecting their AABBs
                double t_l = 0.0, t_r = 0.0;
                {
                    Hit th; BoundingBox tmp(lmin, lmax); if (tmp.intersect(ray, th) && th.hasDistanceAlongRay()) t_l = th.getDistanceAlongRay(); else t_l = std::numeric_limits<double>::infinity();
                }
                {
                    Hit th; BoundingBox tmp(rmin, rmax); if (tmp.intersect(ray, th) && th.hasDistanceAlongRay()) t_r = th.getDistanceAlongRay(); else t_r = std::numeric_limits<double>::infinity();
                }
                if (t_l < t_r) { if (right) st.push({ right, t_r }); if (left) st.push({ left, t_l }); }
                else { if (left) st.push({ left, t_l }); if (right) st.push({ right, t_r }); }
            } else {
                if (left) st.push({ left, 0.0 });
                if (right) st.push({ right, 0.0 });
            }
            continue;
        }

        // Leaf: actual scene object
        Hit h;
        if (node->intersect(ray, h) && h.hasDistanceAlongRay()) {
            const double dist = h.getDistanceAlongRay();
            if (dist < best_dist) { best_dist = dist; best_hit = h; any_hit = true; }
        }
    }
    if (any_hit) out_hit = best_hit;
    return any_hit;
}

double time_trace_all_one_pass(const std::string& json_filepath, const std::string& camera_id) {
	RayTracer rt;
	if (!rt.load_from_json(json_filepath)) return 0.0;
	std::vector<Ray> rays = rt.get_one_ray_per_pixel(camera_id);
	auto t0 = std::chrono::high_resolution_clock::now();
	Hit h;
	for (const Ray& r : rays) {
		(void)rt.one_pass_intersection(r, h);
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dt = t1 - t0;
	return dt.count();
}

double time_trace_all_accelerated_one_pass(const std::string& json_filepath, const std::string& camera_id) {
    RayTracer rt;
    if (!rt.load_from_json(json_filepath)) return 0.0;
    std::vector<Ray> rays = rt.get_one_ray_per_pixel(camera_id);

    // Build object list and BVH
    std::vector<const Mesh*> objects = rt.get_scene_mesh_ptrs();
    std::vector<std::unique_ptr<BoundingBox>> node_store; node_store.reserve(objects.size()*2);
    const Mesh* root = get_BVH_tree(objects, node_store);

    auto t0 = std::chrono::high_resolution_clock::now();
    Hit h;
    for (const Ray& r : rays) {
        (void)accelerated_one_pass_intersection(root, r, h);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    return dt.count();
}
