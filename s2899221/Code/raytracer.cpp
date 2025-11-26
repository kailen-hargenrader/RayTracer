#include "raytracer.h"

#include <fstream>
#include <sstream>
#include <regex>
#include <limits>
#include <chrono>
#include <algorithm>
#include <stack>
#include <cmath>
#include <random>

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
		// Optional camera up vector to preserve roll if present
		double ux=0, uy=0, uz=0;
		if (parse_vec3(cam_block, "up", ux, uy, uz)) {
			Vec3 up(ux, uy, uz);
			cam.updateRotationMatrix(location, direction, up);
		} else {
			cam.updateRotationMatrix(location, direction);
		}
		out[it->str(1)] = cam;
	}
	return !out.empty();
}

static void read_all_meshes(const std::string& content, std::vector<Cube>& cubes, std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres) {
	std::string mesh_block; if (!extract_object_block(content, "MESH", mesh_block)) return;
	std::string cube_block; if (extract_object_block(mesh_block, "Cube", cube_block)) {
		auto v = Cube::read_from_json(cube_block);
		cubes.insert(cubes.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
	}
	std::string plane_block; if (extract_object_block(mesh_block, "Plane", plane_block)) {
		auto v = Plane::read_from_json(plane_block);
		planes.insert(planes.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
	}
	std::string cyl_block; if (extract_object_block(mesh_block, "Cylinder", cyl_block)) {
		auto v = Cylinder::read_from_json(cyl_block);
		cylinders.insert(cylinders.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
	}
	std::string sph_block; if (extract_object_block(mesh_block, "Sphere", sph_block)) {
		auto v = Sphere::read_from_json(sph_block);
		spheres.insert(spheres.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
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
			h.setMesh(m);
			double dist = h.getDistanceAlongRay();
			if (dist < best_dist) { best_dist = dist; best_hit = h; any_hit = true; }
		}
	}
	if (any_hit) { out_hit = best_hit; }
	return any_hit;
}

Pixel RayTracer::shade(const Hit& hit, const Ray& view_ray) const {
	// Return black if no valid surface data
	if (!hit.hasDistanceAlongRay() || !hit.hasSurfaceNormal() || !hit.hasIntersectionPoint()) {
		return Pixel{ 0, 0, 0 };
	}

	// Fetch shading inputs
	const HitVec3 p = hit.getIntersectionPoint();
	HitVec3 n = hit.getSurfaceNormal();

	// Normalize normal defensively
	auto length3 = [](double x, double y, double z) -> double {
		return std::sqrt(x*x + y*y + z*z);
	};
	auto normalize3 = [&](double x, double y, double z) -> HitVec3 {
		const double len = length3(x, y, z);
		if (len <= 1e-15) return HitVec3{ 0.0, 0.0, 1.0 };
		return HitVec3{ x/len, y/len, z/len };
	};
	auto dot3 = [](double ax, double ay, double az, double bx, double by, double bz) -> double {
		return ax*bx + ay*by + az*bz;
	};

	n = normalize3(n.x, n.y, n.z);

	// View direction (from point toward camera)
	const RayVec3 cam_pos = view_ray.getPosition();
	HitVec3 v = normalize3(cam_pos.x - p.x, cam_pos.y - p.y, cam_pos.z - p.z);

	// Material/light parameters
	const double ambient_strength = 0.05;
	const double diffuse_coeff = 1.0;
	const double specular_coeff = 0.5;
	const double shininess = 50.0;

	// Accumulate Blinn-Phong over all point lights (inverse-square falloff)
	double lighting = 0.0;

	if (m_point_lights.empty()) {
		// No lights: renderer will handle producing an all-black image,
		// but keep shading robust and return black here as well.
		//SHOULD NOT GET HERE
		return Pixel{ 0, 0, 0 };
	}
	// Optional metallic adjustment (Principled-style diffuse reduction for metals)
	double metallic_factor = 0.0;
	if (hit.hasMesh()) {
		const Mesh* m = hit.getMesh();
		metallic_factor = std::clamp(m->getMetallic(), 0.0, 1.0);
	}

	for (const auto& kv : m_point_lights) {
		const PointLight& Ls = kv.second;
		const double lx = Ls.location.x - p.x;
		const double ly = Ls.location.y - p.y;
		const double lz = Ls.location.z - p.z;
		const double dist2 = std::max(1e-12, lx*lx + ly*ly + lz*lz);
		const HitVec3 ldir = normalize3(lx, ly, lz);

		// Hard shadow test: cast a ray toward the light, skip if occluded
		{
			const double light_dist = std::sqrt(dist2);
			const double eps = 1e-4;
			// Offset origin along the geometric normal, oriented toward the light, to avoid self-intersection
			const double n_dot_l_tmp = n.x*ldir.x + n.y*ldir.y + n.z*ldir.z;
			const double s = (n_dot_l_tmp >= 0.0) ? 1.0 : -1.0;
			RayVec3 sh_origin{ p.x + s * n.x * eps, p.y + s * n.y * eps, p.z + s * n.z * eps };
			RayVec3 sh_dir{ ldir.x, ldir.y, ldir.z };
			Ray shadow_ray(sh_origin, sh_dir);
			Hit sh;
			bool blocked = false;
			if (one_pass_intersection(shadow_ray, sh) && sh.hasDistanceAlongRay()) {
				const double t = sh.getDistanceAlongRay();
				// Treat as blocked if an occluder lies strictly between the point and the light
				if (t > eps && t < light_dist - eps) {
					blocked = true;
				}
			}
			if (blocked) {
				continue;
			}
		}

		const double n_dot_l = std::max(0.0, dot3(n.x, n.y, n.z, ldir.x, ldir.y, ldir.z));
		// Metals have reduced/zero diffuse in Principled; approximate by (1 - metallic)
		const double diffuse = diffuse_coeff * (1.0 - metallic_factor) * n_dot_l;

		// Blinn-Phong half vector
		const HitVec3 h = normalize3(ldir.x + v.x, ldir.y + v.y, ldir.z + v.z);
		const double n_dot_h = std::max(0.0, dot3(n.x, n.y, n.z, h.x, h.y, h.z));
		// Map roughness to a Blinn-Phong exponent; lower roughness => tighter, stronger highlight
		double roughness = 0.5;
		if (hit.hasMesh()) {
			const Mesh* m = hit.getMesh();
			roughness = std::clamp(m->getRoughness(), 0.001, 1.0);
		}
		const double phong_exp = std::max(2.0, (2.0 / (roughness * roughness)) - 2.0);
		const double specular = (specular_coeff * 1.25) * std::pow(n_dot_h, phong_exp);

		const double attenuation = Ls.radiant_intensity / dist2;
		lighting += attenuation * (diffuse + specular);
	}

	// Add ambient term (reduced for metals) and clamp
	double intensity = (1.0 - metallic_factor) * ambient_strength + lighting;
	if (intensity < 0.0) intensity = 0.0;
	if (intensity > 1.0) intensity = 1.0;

	// Sample albedo from mesh/material
	Pixel base = Pixel{255,255,255};
	if (hit.hasMesh()) {
		const Mesh* m = hit.getMesh();
		base = m->evaluate_albedo(hit);
	}
	auto scale_channel = [&](unsigned char c) -> unsigned char {
		const double cc = intensity * static_cast<double>(c);
		const int ci = static_cast<int>(std::round(cc));
		return static_cast<unsigned char>(std::clamp(ci, 0, 255));
	};
	return Pixel{ scale_channel(base.r), scale_channel(base.g), scale_channel(base.b) };
}

Pixel RayTracer::trace_ray_recursive(const Ray& ray, int depth, int max_scatters, double min_scatter_intensity) const {
	Hit h; (void)one_pass_intersection(ray, h);
	Pixel local = this->shade(h, ray);

	// Stop if maximum depth or no surface
	if (depth >= max_scatters) return local;
	if (!h.hasMesh()) return local;

	// Fresnel-based reflection weight (Schlick), approximating Principled BSDF behavior
	const Mesh* m = h.getMesh();
	const double metallic = m->getMetallic();
	const double roughness = m->getRoughness();
	const double ior = m->getIndexOfRefraction();
	const HitVec3 n = h.getSurfaceNormal();
	const RayVec3 rd = ray.getDirection();
	const double rd_len = std::sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
	const double inv_len = (rd_len > 1e-15) ? (1.0 / rd_len) : 1.0;
	const double vx = -rd.x * inv_len, vy = -rd.y * inv_len, vz = -rd.z * inv_len;
	const double cosTheta = std::max(0.0, n.x*vx + n.y*vy + n.z*vz);
	// F0: dielectrics from IOR, metals from base color intensity
	double F0 = 0.04;
	if (metallic <= 0.0) {
		const double r0 = (ior - 1.0) / (ior + 1.0);
		F0 = std::clamp(r0 * r0, 0.0, 1.0);
	} else {
		Pixel base = m->evaluate_albedo(h);
		const double avg = (static_cast<double>(base.r) + static_cast<double>(base.g) + static_cast<double>(base.b)) / (3.0 * 255.0);
		F0 = std::clamp(avg, 0.0, 1.0);
	}
	double Fr = F0 + (1.0 - F0) * std::pow(1.0 - cosTheta, 5.0);
	// Keep Fresnel magnitude; roughness affects lobe shape, not energy here

	// Skip very weak contributions
	if (Fr <= 0.0 || Fr < min_scatter_intensity) {
		return local;
	}

	// Build secondary ray with small offset to avoid self-intersection
	Ray secondary_ray = ray; // will overwrite below
	{
		const HitVec3 ip = h.getIntersectionPoint();
		RayVec3 sec_origin{ ip.x, ip.y, ip.z };
		RayVec3 sec_dir = ray.getDirection();
		if (h.hasReflectedDirection()) {
			const HitVec3 r = h.getReflectedDirection();
			sec_dir = RayVec3{ r.x, r.y, r.z };
		} else {
			// r = d - 2 (d·n) n
			const double d_dot_n = rd.x*n.x + rd.y*n.y + rd.z*n.z;
			sec_dir = RayVec3{ rd.x - 2.0*d_dot_n*n.x, rd.y - 2.0*d_dot_n*n.y, rd.z - 2.0*d_dot_n*n.z };
		}
		const double eps = 1e-6;
		sec_origin.x += sec_dir.x * eps;
		sec_origin.y += sec_dir.y * eps;
		sec_origin.z += sec_dir.z * eps;
		secondary_ray = Ray(sec_origin, sec_dir);
	}

	// Recurse; misses contribute black naturally via shade()
	const Pixel sec_color = trace_ray_recursive(secondary_ray, depth + 1, max_scatters, min_scatter_intensity);

	// For metals, tint reflections by base color; for dielectrics keep reflections achromatic
	Pixel sec_tinted = sec_color;
	if (metallic > 0.0) {
		Pixel base = m->evaluate_albedo(h);
		auto mul = [](unsigned char a, unsigned char b) -> unsigned char {
			const double aa = static_cast<double>(a) / 255.0;
			const double bb = static_cast<double>(b);
			const int ci = static_cast<int>(std::round(aa * bb));
			return static_cast<unsigned char>(std::clamp(ci, 0, 255));
		};
		sec_tinted.r = mul(base.r, sec_color.r);
		sec_tinted.g = mul(base.g, sec_color.g);
		sec_tinted.b = mul(base.b, sec_color.b);
	}

	// Blend local and reflection by Fresnel
	auto blend_channel = [&](unsigned char a, unsigned char b) -> unsigned char {
		const double aa = static_cast<double>(a);
		const double bb = static_cast<double>(b);
		const double cc = (1.0 - Fr) * aa + Fr * bb;
		const int ci = static_cast<int>(std::round(cc));
		return static_cast<unsigned char>(std::clamp(ci, 0, 255));
	};

	Pixel out;
	out.r = blend_channel(local.r, sec_tinted.r);
	out.g = blend_channel(local.g, sec_tinted.g);
	out.b = blend_channel(local.b, sec_tinted.b);
	return out;
}

bool RayTracer::render_unaccelerated_ppm(const std::string& camera_id, const std::string& output_filepath, int samples_per_pixel, int max_scatters, double min_scatter_intensity) const {
	auto it = m_cameras.find(camera_id);
	if (it == m_cameras.end()) return false;
	const Camera& cam = it->second;

	int width = 0, height = 0;
	cam.getImageResolution(width, height);
	if (width <= 0 || height <= 0) return false;
	if (samples_per_pixel < 1) samples_per_pixel = 1;

	Image img(width, height);
	img.setMaxValue(255);

	// If no lights in the scene, produce an all-black image immediately
	if (m_point_lights.empty()) {
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				img.setPixel(x, y, 0, 0, 0);
			}
		}
		img.write(output_filepath);
		return true;
	}

	// Random jitter in pixel space, bounded to half a pixel in each axis
	std::mt19937 rng(static_cast<unsigned int>(std::random_device{}()));
	std::uniform_real_distribution<double> jitter_dist(-0.5, 0.5);

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			double accum_r = 0.0;
			double accum_g = 0.0;
			double accum_b = 0.0;

			for (int s = 0; s < samples_per_pixel; ++s) {
				const double jx = (s == 0) ? 0.0 : jitter_dist(rng);
				const double jy = (s == 0) ? 0.0 : jitter_dist(rng);

				Vec3 ro, rd;
				const double px = static_cast<double>(x) + 0.5 + jx;
				const double py = static_cast<double>(y) + 0.5 + jy;
				cam.pixelToRay(px, py, ro, rd);

				Ray ray(RayVec3{ ro.x, ro.y, ro.z }, RayVec3{ rd.x, rd.y, rd.z });
				Pixel spx = trace_ray_recursive(ray, 0, max_scatters, min_scatter_intensity);

				accum_r += static_cast<double>(spx.r);
				accum_g += static_cast<double>(spx.g);
				accum_b += static_cast<double>(spx.b);
			}

			const int rr = static_cast<int>(std::round(accum_r / static_cast<double>(samples_per_pixel)));
			const int gg = static_cast<int>(std::round(accum_g / static_cast<double>(samples_per_pixel)));
			const int bb = static_cast<int>(std::round(accum_b / static_cast<double>(samples_per_pixel)));
			// Gamma encode to sRGB for display (Blender uses view transform; approximate with gamma 2.2)
			auto gamma_encode = [](int c_lin) -> int {
				double v = std::clamp(static_cast<double>(c_lin), 0.0, 255.0) / 255.0;
				double v_srgb = std::pow(v, 1.0 / 2.2);
				return static_cast<int>(std::round(v_srgb * 255.0));
			};
			const int r_out = gamma_encode(rr);
			const int g_out = gamma_encode(gg);
			const int b_out = gamma_encode(bb);
			img.setPixel(x, y, r_out, g_out, b_out);
		}
	}

	img.write(output_filepath);
	return true;
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
			h.setMesh(node);
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
