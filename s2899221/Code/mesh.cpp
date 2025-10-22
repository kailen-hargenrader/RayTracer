#include "mesh.h"
#include "utils.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <cmath>
#include <limits>

// JSON-like parsing helpers are implemented in utils.{h,cpp}

// ---------- Cube ----------
Cube::Cube() : m_translation{0,0,0}, m_rotation{0,0,0}, m_scale(1.0) {}
Cube::Cube(const Float3& translation, const EulerAngles& rotation, double scale)
    : m_translation(translation), m_rotation(rotation), m_scale(scale) {}

std::vector<Cube> Cube::read_from_json(const std::string& class_block) {
    // class_block expected to contain entries keyed by id: { "1": { ... }, "2": { ... }, ... }
    std::vector<Cube> result;

    // Find each object id by regex on quoted number keys
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        // Extract the object block following this id
        std::string sub;
        // Build a temporary string starting from this position to ensure extract works
        std::string tail = class_block.substr(key_pos);
        if (!util_json::extract_object_block(tail, it->str(1), sub)) continue;

        Float3 tr{0,0,0}; EulerAngles rot{0,0,0}; double scale = 1.0;
        {
            double x=0,y=0,z=0; if (util_json::parse_vec3(sub, "translation", x, y, z)) tr = Float3{ x, y, z };
        }
        // rotation fields named roll, pitch, yaw
        double r=0,p=0,y=0; util_json::parse_number(sub, "roll", r); util_json::parse_number(sub, "pitch", p); util_json::parse_number(sub, "yaw", y);
        rot = EulerAngles{ r, p, y };
        util_json::parse_number(sub, "scale", scale);

        result.emplace_back(tr, rot, scale);
    }
    return result;
}

void Cube::write_to_console(std::ostream& out) const {
    out << "Cube(translation=[" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z
        << "], rotation(rpy)=[" << m_rotation.roll << ", " << m_rotation.pitch << ", " << m_rotation.yaw
        << "], scale=" << m_scale << ")\n";
}

static inline void build_rotation_matrix_rpy(double roll, double pitch, double yaw, double R[3][3]) {
    const double cr = std::cos(roll), sr = std::sin(roll);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double cy = std::cos(yaw), sy = std::sin(yaw);

    // R = Rz(yaw) * Ry(pitch) * Rx(roll)
    const double Rz[3][3] = { { cy, -sy, 0 }, { sy, cy, 0 }, { 0, 0, 1 } };
    const double Ry[3][3] = { { cp, 0, sp }, { 0, 1, 0 }, { -sp, 0, cp } };
    const double Rx[3][3] = { { 1, 0, 0 }, { 0, cr, -sr }, { 0, sr, cr } };

    double A[3][3];
    // A = Ry * Rx
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j] = Ry[i][0]*Rx[0][j] + Ry[i][1]*Rx[1][j] + Ry[i][2]*Rx[2][j];
        }
    }
    // R = Rz * A
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            R[i][j] = Rz[i][0]*A[0][j] + Rz[i][1]*A[1][j] + Rz[i][2]*A[2][j];
        }
    }
}

static inline void mul_mat_vec(const double M[3][3], const double v[3], double out[3]) {
    out[0] = M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2];
    out[1] = M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2];
    out[2] = M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2];
}

static inline void transpose3(const double M[3][3], double MT[3][3]) {
    MT[0][0] = M[0][0]; MT[0][1] = M[1][0]; MT[0][2] = M[2][0];
    MT[1][0] = M[0][1]; MT[1][1] = M[1][1]; MT[1][2] = M[2][1];
    MT[2][0] = M[0][2]; MT[2][1] = M[1][2]; MT[2][2] = M[2][2];
}

bool Cube::intersect(const Ray& ray, Hit& hit) const {
    // Build rotation matrix for local->world; we need its transpose for world->local
    double R[3][3];
    build_rotation_matrix_rpy(m_rotation.roll, m_rotation.pitch, m_rotation.yaw, R);
    double RT[3][3];
    transpose3(R, RT);

    // World ray
    const RayVec3 ro_w = ray.getPosition();
    const RayVec3 rd_w = ray.getDirection();

    // Translate to object frame and apply inverse rotation and inverse uniform scale
    const double to_obj[3] = { ro_w.x - m_translation.x, ro_w.y - m_translation.y, ro_w.z - m_translation.z };
    double ro_r[3]; mul_mat_vec(RT, to_obj, ro_r);
    double rd_in[3] = { rd_w.x, rd_w.y, rd_w.z };
    double rd_r[3]; mul_mat_vec(RT, rd_in, rd_r);

    const double invS = (m_scale != 0.0) ? (1.0 / m_scale) : 0.0;
    double ro[3] = { ro_r[0] * invS, ro_r[1] * invS, ro_r[2] * invS };
    double rd[3] = { rd_r[0] * invS, rd_r[1] * invS, rd_r[2] * invS };

    // Slab intersection with AABB [-1,1]^3
    const double eps = 1e-12;
    double tmin = -std::numeric_limits<double>::infinity();
    double tmax =  std::numeric_limits<double>::infinity();
    double n_enter[3] = {0,0,0};
    double n_exit[3]  = {0,0,0};

    for (int axis = 0; axis < 3; ++axis) {
        const double o = ro[axis];
        const double d = rd[axis];
        const double low = -1.0, high = 1.0;

        if (std::abs(d) < eps) {
            if (o < low || o > high) return false; // parallel and outside
            // parallel and inside; no update to tmin/tmax
            continue;
        }

        double t1 = (low - o) / d;
        double t2 = (high - o) / d;

        double n1[3] = {0,0,0}; n1[axis] = -1.0; // normal for low plane
        double n2[3] = {0,0,0}; n2[axis] =  1.0; // normal for high plane

        if (t1 > t2) { std::swap(t1, t2); std::swap(n1, n2); }

        if (t1 > tmin) { tmin = t1; n_enter[0] = n1[0]; n_enter[1] = n1[1]; n_enter[2] = n1[2]; }
        if (t2 < tmax) { tmax = t2; n_exit[0]  = n2[0]; n_exit[1]  = n2[1]; n_exit[2]  = n2[2]; }

        if (tmin > tmax) return false;
    }

    if (tmax < 0.0) return false; // box is behind the ray origin

    const double t_hit = (tmin >= 0.0) ? tmin : tmax; // if inside, use exit
    const bool used_enter = (tmin >= 0.0);
    const double* n_local = used_enter ? n_enter : n_exit;

    // Intersection point in local cube space
    double p_local[3] = { ro[0] + t_hit * rd[0], ro[1] + t_hit * rd[1], ro[2] + t_hit * rd[2] };

    // Transform normal back to world (rotation only)
    double n_world[3]; mul_mat_vec(R, n_local, n_world);
    // Normalize normal to be unit length
    const double nlen = std::sqrt(n_world[0]*n_world[0] + n_world[1]*n_world[1] + n_world[2]*n_world[2]);
    double n_unit[3] = { n_world[0]/nlen, n_world[1]/nlen, n_world[2]/nlen };

    // Transform intersection point back to world: p_w = T + R * (p_local * scale)
    double p_scaled[3] = { p_local[0] * m_scale, p_local[1] * m_scale, p_local[2] * m_scale };
    double p_world[3]; mul_mat_vec(R, p_scaled, p_world);
    p_world[0] += m_translation.x; p_world[1] += m_translation.y; p_world[2] += m_translation.z;

    // Reflected direction in world space: r = d - 2 (d·n) n (n must be unit)
    const double d_world[3] = { rd_w.x, rd_w.y, rd_w.z };
    const double d_dot_n = d_world[0]*n_unit[0] + d_world[1]*n_unit[1] + d_world[2]*n_unit[2];
    double r_world[3] = { d_world[0] - 2.0 * d_dot_n * n_unit[0],
                          d_world[1] - 2.0 * d_dot_n * n_unit[1],
                          d_world[2] - 2.0 * d_dot_n * n_unit[2] };

    // Populate Hit (world coordinates)
    HitVec3 ip{ p_world[0], p_world[1], p_world[2] };
    HitVec3 nrm{ n_unit[0], n_unit[1], n_unit[2] };
    HitVec3 refl{ r_world[0], r_world[1], r_world[2] };

    hit.setIntersectionPoint(ip);
    hit.setSurfaceNormal(nrm);
    hit.setReflectedDirection(refl);
    // Distance along original world ray
    const double dx = p_world[0] - ro_w.x;
    const double dy = p_world[1] - ro_w.y;
    const double dz = p_world[2] - ro_w.z;
    const double dist_world = std::sqrt(dx*dx + dy*dy + dz*dz);
    hit.setDistanceAlongRay(dist_world);
    return true;
}

// ---------- Cylinder ----------
Cylinder::Cylinder() : m_translation{0,0,0}, m_rotation{0,0,0}, m_scale(1.0), m_radius(1.0), m_length(1.0) {}
Cylinder::Cylinder(const Float3& translation, const EulerAngles& rotation, double scale, double radius, double length)
    : m_translation(translation), m_rotation(rotation), m_scale(scale), m_radius(radius), m_length(length) {}

std::vector<Cylinder> Cylinder::read_from_json(const std::string& class_block) {
    std::vector<Cylinder> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!util_json::extract_object_block(tail, it->str(1), sub)) continue;

        Float3 tr{0,0,0}; EulerAngles rot{0,0,0}; double scale = 1.0; double radius = 1.0; double length = 1.0;
        { double x=0,y=0,z=0; if (util_json::parse_vec3(sub, "translation", x, y, z)) tr = Float3{ x, y, z }; }
        double r=0,p=0,y=0; util_json::parse_number(sub, "roll", r); util_json::parse_number(sub, "pitch", p); util_json::parse_number(sub, "yaw", y);
        rot = EulerAngles{ r, p, y };
        util_json::parse_number(sub, "scale", scale);
        util_json::parse_number(sub, "radius", radius);
        util_json::parse_number(sub, "length", length);

        result.emplace_back(tr, rot, scale, radius, length);
    }
    return result;
}

void Cylinder::write_to_console(std::ostream& out) const {
    out << "Cylinder(translation=[" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z
        << "], rotation(rpy)=[" << m_rotation.roll << ", " << m_rotation.pitch << ", " << m_rotation.yaw
        << "], scale=" << m_scale << ", radius=" << m_radius << ", length=" << m_length << ")\n";
}

// ---------- Sphere ----------
Sphere::Sphere() : m_location{0,0,0}, m_radius(1.0) {}
Sphere::Sphere(const Float3& location, double radius) : m_location(location), m_radius(radius) {}

std::vector<Sphere> Sphere::read_from_json(const std::string& class_block) {
    std::vector<Sphere> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!util_json::extract_object_block(tail, it->str(1), sub)) continue;

        Float3 loc{0,0,0}; double radius = 1.0;
        { double x=0,y=0,z=0; if (util_json::parse_vec3(sub, "location", x, y, z)) loc = Float3{ x, y, z }; }
        util_json::parse_number(sub, "radius", radius);
        result.emplace_back(loc, radius);
    }
    return result;
}

void Sphere::write_to_console(std::ostream& out) const {
    out << "Sphere(location=[" << m_location.x << ", " << m_location.y << ", " << m_location.z
        << "], radius=" << m_radius << ")\n";
}

// ---------- Plane ----------
Plane::Plane() : m_corners() {}
Plane::Plane(const std::vector<Float3>& corners) : m_corners(corners) {}

std::vector<Plane> Plane::read_from_json(const std::string& class_block) {
    std::vector<Plane> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!util_json::extract_object_block(tail, it->str(1), sub)) continue;

        // corners is an array of 4 objects {x,y,z}
        std::string corners_block;
        if (!util_json::extract_object_block(sub, "corners", corners_block)) {
            continue;
        }
        // corners_block here is content between '[' and ']' thanks to support in extract_object_block
        // Parse up to 4 corner objects in order of appearance
        std::vector<Float3> corners;
        size_t search_offset = 0;
        while (true) {
            const size_t brace_pos = corners_block.find('{', search_offset);
            if (brace_pos == std::string::npos) break;

            // Manual brace match to find end of this corner object
            int depth = 0; bool in_str = false; bool esc = false; size_t end_pos = std::string::npos;
            for (size_t j = brace_pos; j < corners_block.size(); ++j) {
                char c = corners_block[j];
                if (in_str) { if (esc) esc = false; else if (c == '\\') esc = true; else if (c == '"') in_str = false; continue; }
                if (c == '"') { in_str = true; continue; }
                if (c == '{') { ++depth; }
                else if (c == '}') { --depth; if (depth == 0) { end_pos = j; break; } }
            }
            if (end_pos == std::string::npos) break;

            const std::string one_corner = corners_block.substr(brace_pos + 1, end_pos - brace_pos - 1);
            search_offset = end_pos + 1;

            double x=0,y=0,z=0; util_json::parse_number(one_corner, "x", x); util_json::parse_number(one_corner, "y", y); util_json::parse_number(one_corner, "z", z);
            corners.push_back(Float3{ x, y, z });
            if (corners.size() == 4) break;
        }
        if (!corners.empty()) {
            result.emplace_back(corners);
        }
    }
    return result;
}

void Plane::write_to_console(std::ostream& out) const {
    out << "Plane(corners=[";
    for (size_t i = 0; i < m_corners.size(); ++i) {
        const auto& c = m_corners[i];
        out << "[" << c.x << ", " << c.y << ", " << c.z << "]";
        if (i + 1 < m_corners.size()) out << ", ";
    }
    out << "])\n";
}


