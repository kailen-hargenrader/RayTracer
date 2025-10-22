#include "mesh.h"
#include "utils.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <cmath>
#include <limits>
#include <algorithm>

// JSON-like parsing helpers are implemented in utils.{h,cpp}

// ---------- Cube ----------
Cube::Cube() : m_translation{0,0,0}, m_rotation{0,0,0}, m_scale{1.0, 1.0, 1.0} {}
Cube::Cube(const Float3& translation, const EulerAngles& rotation, const Float3& scale)
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

        Float3 tr{0,0,0}; EulerAngles rot{0,0,0}; Float3 scale{1.0, 1.0, 1.0};
        {
            double x=0,y=0,z=0; if (util_json::parse_vec3(sub, "translation", x, y, z)) tr = Float3{ x, y, z };
        }
        // rotation fields named roll, pitch, yaw
        double r=0,p=0,y=0; util_json::parse_number(sub, "roll", r); util_json::parse_number(sub, "pitch", p); util_json::parse_number(sub, "yaw", y);
        rot = EulerAngles{ r, p, y };
        {
            double sx=1, sy=1, sz=1; if (util_json::parse_vec3(sub, "scale", sx, sy, sz)) scale = Float3{ sx, sy, sz };
        }

        result.emplace_back(tr, rot, scale);
    }
    return result;
}

void Cube::write_to_console(std::ostream& out) const {
    out << "Cube(translation=[" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z
        << "], rotation(rpy)=[" << m_rotation.roll << ", " << m_rotation.pitch << ", " << m_rotation.yaw
        << "], scale=[" << m_scale.x << ", " << m_scale.y << ", " << m_scale.z << "])\n";
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

    // Translate to object frame and apply inverse rotation and inverse non-uniform scale
    const double to_obj[3] = { ro_w.x - m_translation.x, ro_w.y - m_translation.y, ro_w.z - m_translation.z };
    double ro_r[3]; mul_mat_vec(RT, to_obj, ro_r);
    double rd_in[3] = { rd_w.x, rd_w.y, rd_w.z };
    double rd_r[3]; mul_mat_vec(RT, rd_in, rd_r);

    const double invSx = (m_scale.x != 0.0) ? (1.0 / m_scale.x) : 0.0;
    const double invSy = (m_scale.y != 0.0) ? (1.0 / m_scale.y) : 0.0;
    const double invSz = (m_scale.z != 0.0) ? (1.0 / m_scale.z) : 0.0;
    double ro[3] = { ro_r[0] * invSx, ro_r[1] * invSy, ro_r[2] * invSz };
    double rd[3] = { rd_r[0] * invSx, rd_r[1] * invSy, rd_r[2] * invSz };

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

    // Transform normal back to world: n_world ∝ R * S^{-1} * n_local
    double n_scaled[3] = { n_local[0] * invSx, n_local[1] * invSy, n_local[2] * invSz };
    double n_world[3]; mul_mat_vec(R, n_scaled, n_world);
    // Normalize normal to be unit length
    const double nlen = std::sqrt(n_world[0]*n_world[0] + n_world[1]*n_world[1] + n_world[2]*n_world[2]);
    double n_unit[3] = { n_world[0]/nlen, n_world[1]/nlen, n_world[2]/nlen };

    // Transform intersection point back to world: p_w = T + R * (S * p_local)
    double p_scaled[3] = { p_local[0] * m_scale.x, p_local[1] * m_scale.y, p_local[2] * m_scale.z };
    double p_world[3]; mul_mat_vec(R, p_scaled, p_world);
    p_world[0] += m_translation.x; p_world[1] += m_translation.y; p_world[2] += m_translation.z;

    // Reflected direction in world space: r = d - 2 (d·n) n (ensure unit length)
    const double d_world[3] = { rd_w.x, rd_w.y, rd_w.z };
    const double d_dot_n = d_world[0]*n_unit[0] + d_world[1]*n_unit[1] + d_world[2]*n_unit[2];
    double r_world_tmp[3] = { d_world[0] - 2.0 * d_dot_n * n_unit[0],
                              d_world[1] - 2.0 * d_dot_n * n_unit[1],
                              d_world[2] - 2.0 * d_dot_n * n_unit[2] };
    double r_len = std::sqrt(r_world_tmp[0]*r_world_tmp[0] + r_world_tmp[1]*r_world_tmp[1] + r_world_tmp[2]*r_world_tmp[2]);
    if (r_len < 1e-15) r_len = 1.0; // avoid divide by zero
    double r_world[3] = { r_world_tmp[0]/r_len, r_world_tmp[1]/r_len, r_world_tmp[2]/r_len };

    // Populate Hit (world coordinates)
    HitVec3 ip{ p_world[0], p_world[1], p_world[2] };
    HitVec3 nrm{ n_unit[0], n_unit[1], n_unit[2] };
    HitVec3 refl{ r_world[0], r_world[1], r_world[2] };

    hit.setIntersectionPoint(ip);
    hit.setSurfaceNormal(nrm);
    hit.setReflectedDirection(refl);
    // Distance along original world ray (projection along the ray direction)
    const double dx = p_world[0] - ro_w.x;
    const double dy = p_world[1] - ro_w.y;
    const double dz = p_world[2] - ro_w.z;
    const double d_len = std::sqrt(d_world[0]*d_world[0] + d_world[1]*d_world[1] + d_world[2]*d_world[2]);
    const double dist_world = (d_len > 1e-15) ? (dx*d_world[0] + dy*d_world[1] + dz*d_world[2]) / d_len
                                              : std::sqrt(dx*dx + dy*dy + dz*dz);
    hit.setDistanceAlongRay(dist_world);
    return true;
}

static inline void aabb_expand(Float3& bmin, Float3& bmax, const double p[3]) {
    if (p[0] < bmin.x) bmin.x = p[0]; if (p[1] < bmin.y) bmin.y = p[1]; if (p[2] < bmin.z) bmin.z = p[2];
    if (p[0] > bmax.x) bmax.x = p[0]; if (p[1] > bmax.y) bmax.y = p[1]; if (p[2] > bmax.z) bmax.z = p[2];
}

// Compute tight world AABB for transformed cube by transforming its 8 corners
void Cube::compute_aabb(Float3& outMin, Float3& outMax) const {
    double R[3][3];
    build_rotation_matrix_rpy(m_rotation.roll, m_rotation.pitch, m_rotation.yaw, R);
    auto transform_point = [&](double lx, double ly, double lz, double out[3]) {
        double p_scaled[3] = { lx * m_scale.x, ly * m_scale.y, lz * m_scale.z };
        double p_rot[3]; mul_mat_vec(R, p_scaled, p_rot);
        out[0] = p_rot[0] + m_translation.x;
        out[1] = p_rot[1] + m_translation.y;
        out[2] = p_rot[2] + m_translation.z;
    };
    outMin = { +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity() };
    outMax = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
    for (int dx = -1; dx <= 1; dx += 2) {
        for (int dy = -1; dy <= 1; dy += 2) {
            for (int dz = -1; dz <= 1; dz += 2) {
                double pw[3];
                transform_point(static_cast<double>(dx), static_cast<double>(dy), static_cast<double>(dz), pw);
                aabb_expand(outMin, outMax, pw);
            }
        }
    }
}

// ---------- Cylinder ----------
Cylinder::Cylinder() : m_translation{0,0,0}, m_rotation{0,0,0}, m_scale{1.0,1.0,1.0} {}
Cylinder::Cylinder(const Float3& translation, const EulerAngles& rotation, const Float3& scale)
    : m_translation(translation), m_rotation(rotation), m_scale(scale) {}

std::vector<Cylinder> Cylinder::read_from_json(const std::string& class_block) {
    std::vector<Cylinder> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!util_json::extract_object_block(tail, it->str(1), sub)) continue;

        Float3 tr{0,0,0}; EulerAngles rot{0,0,0}; Float3 scale{1.0,1.0,1.0};
        { double x=0,y=0,z=0; if (util_json::parse_vec3(sub, "translation", x, y, z)) tr = Float3{ x, y, z }; }
        double rr=0,pp=0,yy=0; util_json::parse_number(sub, "roll", rr); util_json::parse_number(sub, "pitch", pp); util_json::parse_number(sub, "yaw", yy);
        rot = EulerAngles{ rr, pp, yy };
        { double sx=1,sy=1,sz=1; if (util_json::parse_vec3(sub, "scale", sx, sy, sz)) scale = Float3{ sx, sy, sz }; }
        result.emplace_back(tr, rot, scale);
    }
    return result;
}

void Cylinder::write_to_console(std::ostream& out) const {
    out << "Cylinder(translation=[" << m_translation.x << ", " << m_translation.y << ", " << m_translation.z
        << "], rotation(rpy)=[" << m_rotation.roll << ", " << m_rotation.pitch << ", " << m_rotation.yaw
        << "], scale=[" << m_scale.x << ", " << m_scale.y << ", " << m_scale.z << "])\n";
}

bool Cylinder::intersect(const Ray& ray, Hit& hit) const {
    // Build rotation matrix
    double R[3][3];
    build_rotation_matrix_rpy(m_rotation.roll, m_rotation.pitch, m_rotation.yaw, R);
    double RT[3][3];
    transpose3(R, RT);

    // World ray
    const RayVec3 ro_w = ray.getPosition();
    const RayVec3 rd_w = ray.getDirection();

    // Transform to local coordinates (non-uniform scale)
    const double to_obj[3] = { ro_w.x - m_translation.x, ro_w.y - m_translation.y, ro_w.z - m_translation.z };
    double ro_r[3]; mul_mat_vec(RT, to_obj, ro_r);
    double rd_in[3] = { rd_w.x, rd_w.y, rd_w.z };
    double rd_r[3]; mul_mat_vec(RT, rd_in, rd_r);

    const double invSx = (m_scale.x != 0.0) ? (1.0 / m_scale.x) : 0.0;
    const double invSy = (m_scale.y != 0.0) ? (1.0 / m_scale.y) : 0.0;
    const double invSz = (m_scale.z != 0.0) ? (1.0 / m_scale.z) : 0.0;
    double ro[3] = { ro_r[0] * invSx, ro_r[1] * invSy, ro_r[2] * invSz };
    double rd[3] = { rd_r[0] * invSx, rd_r[1] * invSy, rd_r[2] * invSz };

    // Canonical cylinder: radius 1 in x-y, z in [-1,1]. World scale stretches it.
    const double r = 1.0;
    const double halfL = 1.0;
    const double eps = 1e-12;

    double best_t = std::numeric_limits<double>::infinity();
    double best_p[3] = {0,0,0};
    double best_n_local[3] = {0,0,0};
    bool hit_found = false;

    // Side intersection: x^2 + y^2 = r^2
    const double A = rd[0]*rd[0] + rd[1]*rd[1];
    const double B = 2.0 * (ro[0]*rd[0] + ro[1]*rd[1]);
    const double C = ro[0]*ro[0] + ro[1]*ro[1] - r*r;
    if (A > eps) {
        const double disc = B*B - 4.0*A*C;
        if (disc >= 0.0) {
            const double sqrt_disc = std::sqrt(disc);
            const double t1 = (-B - sqrt_disc) / (2.0*A);
            const double t2 = (-B + sqrt_disc) / (2.0*A);
            auto try_t = [&](double t) {
                if (t < 0.0) return;
                const double z = ro[2] + t * rd[2];
                if (z < -halfL - eps || z > halfL + eps) return;
                if (t < best_t) {
                    best_t = t;
                    best_p[0] = ro[0] + t*rd[0];
                    best_p[1] = ro[1] + t*rd[1];
                    best_p[2] = z;
                    // normal outward on side
                    const double len_xy = std::sqrt(best_p[0]*best_p[0] + best_p[1]*best_p[1]);
                    if (len_xy > eps) {
                        best_n_local[0] = best_p[0]/len_xy;
                        best_n_local[1] = best_p[1]/len_xy;
                        best_n_local[2] = 0.0;
                    } else {
                        best_n_local[0] = best_n_local[1] = 0.0; best_n_local[2] = 1.0;
                    }
                    hit_found = true;
                }
            };
            try_t(t1);
            try_t(t2);
        }
    }

    // Caps at z = ±halfL
    auto try_cap = [&](double z_cap, double nz) {
        if (std::abs(rd[2]) < eps) return; // parallel to cap plane
        const double t = (z_cap - ro[2]) / rd[2];
        if (t < 0.0) return;
        const double x = ro[0] + t * rd[0];
        const double y = ro[1] + t * rd[1];
        if (x*x + y*y <= r*r + 1e-9) {
            if (t < best_t) {
                best_t = t;
                best_p[0] = x; best_p[1] = y; best_p[2] = z_cap;
                best_n_local[0] = 0.0; best_n_local[1] = 0.0; best_n_local[2] = nz;
                hit_found = true;
            }
        }
    };
    try_cap(+halfL, +1.0);
    try_cap(-halfL, -1.0);

    if (!hit_found) return false;

    // Transform normal back to world: n_world ∝ R * S^{-1} * n_local
    double n_scaled[3] = { best_n_local[0] * invSx, best_n_local[1] * invSy, best_n_local[2] * invSz };
    double n_world_raw[3]; mul_mat_vec(R, n_scaled, n_world_raw);
    const double nlen = std::sqrt(n_world_raw[0]*n_world_raw[0] + n_world_raw[1]*n_world_raw[1] + n_world_raw[2]*n_world_raw[2]);
    double n_unit[3] = { n_world_raw[0]/nlen, n_world_raw[1]/nlen, n_world_raw[2]/nlen };

    double p_scaled[3] = { best_p[0] * m_scale.x, best_p[1] * m_scale.y, best_p[2] * m_scale.z };
    double p_world_rel[3]; mul_mat_vec(R, p_scaled, p_world_rel);
    double p_world[3] = { p_world_rel[0] + m_translation.x, p_world_rel[1] + m_translation.y, p_world_rel[2] + m_translation.z };

    // Reflected direction in world space
    const double d_world[3] = { rd_w.x, rd_w.y, rd_w.z };
    const double d_dot_n = d_world[0]*n_unit[0] + d_world[1]*n_unit[1] + d_world[2]*n_unit[2];
    double r_world_tmp[3] = { d_world[0] - 2.0 * d_dot_n * n_unit[0],
                              d_world[1] - 2.0 * d_dot_n * n_unit[1],
                              d_world[2] - 2.0 * d_dot_n * n_unit[2] };
    double r_len = std::sqrt(r_world_tmp[0]*r_world_tmp[0] + r_world_tmp[1]*r_world_tmp[1] + r_world_tmp[2]*r_world_tmp[2]);
    if (r_len < 1e-15) r_len = 1.0;
    double r_world[3] = { r_world_tmp[0]/r_len, r_world_tmp[1]/r_len, r_world_tmp[2]/r_len };

    HitVec3 ip{ p_world[0], p_world[1], p_world[2] };
    HitVec3 nrm{ n_unit[0], n_unit[1], n_unit[2] };
    HitVec3 refl{ r_world[0], r_world[1], r_world[2] };

    hit.setIntersectionPoint(ip);
    hit.setSurfaceNormal(nrm);
    hit.setReflectedDirection(refl);

    const double dx = p_world[0] - ro_w.x;
    const double dy = p_world[1] - ro_w.y;
    const double dz = p_world[2] - ro_w.z;
    const double d_len = std::sqrt(d_world[0]*d_world[0] + d_world[1]*d_world[1] + d_world[2]*d_world[2]);
    const double dist_world = (d_len > 1e-15) ? (dx*d_world[0] + dy*d_world[1] + dz*d_world[2]) / d_len
                                              : std::sqrt(dx*dx + dy*dy + dz*dz);
    hit.setDistanceAlongRay(dist_world);
    return true;
}

// Conservative AABB via transforming the 8 corners of local cube [-1,1]^3
void Cylinder::compute_aabb(Float3& outMin, Float3& outMax) const {
    double R[3][3];
    build_rotation_matrix_rpy(m_rotation.roll, m_rotation.pitch, m_rotation.yaw, R);
    auto transform_point = [&](double lx, double ly, double lz, double out[3]) {
        double p_scaled[3] = { lx * m_scale.x, ly * m_scale.y, lz * m_scale.z };
        double p_rot[3]; mul_mat_vec(R, p_scaled, p_rot);
        out[0] = p_rot[0] + m_translation.x;
        out[1] = p_rot[1] + m_translation.y;
        out[2] = p_rot[2] + m_translation.z;
    };
    outMin = { +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity() };
    outMax = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
    for (int dx = -1; dx <= 1; dx += 2) {
        for (int dy = -1; dy <= 1; dy += 2) {
            for (int dz = -1; dz <= 1; dz += 2) {
                double pw[3];
                transform_point(static_cast<double>(dx), static_cast<double>(dy), static_cast<double>(dz), pw);
                aabb_expand(outMin, outMax, pw);
            }
        }
    }
}

// ---------- Sphere ----------
Sphere::Sphere() : m_location{0,0,0}, m_scale{1.0,1.0,1.0} {}
Sphere::Sphere(const Float3& location, const Float3& scale) : m_location(location), m_scale(scale) {}

std::vector<Sphere> Sphere::read_from_json(const std::string& class_block) {
    std::vector<Sphere> result;
    std::regex id_rx("\\\"([0-9]+)\\\"\\s*:\\s*\\{");
    auto begin = std::sregex_iterator(class_block.begin(), class_block.end(), id_rx);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        const size_t key_pos = static_cast<size_t>(it->position());
        std::string tail = class_block.substr(key_pos);
        std::string sub; if (!util_json::extract_object_block(tail, it->str(1), sub)) continue;

        Float3 loc{0,0,0}; Float3 scale{1.0,1.0,1.0};
        { double x=0,y=0,z=0; if (util_json::parse_vec3(sub, "location", x, y, z)) loc = Float3{ x, y, z }; }
        // Prefer per-axis scale if provided; else radius for uniform sphere
        double sx=1, sy=1, sz=1; if (util_json::parse_vec3(sub, "scale", sx, sy, sz)) { scale = Float3{ sx, sy, sz }; }
        result.emplace_back(loc, scale);
    }
    return result;
}

void Sphere::write_to_console(std::ostream& out) const {
    out << "Sphere(location=[" << m_location.x << ", " << m_location.y << ", " << m_location.z
        << "], scale=[" << m_scale.x << ", " << m_scale.y << ", " << m_scale.z << "])\n";
}

bool Sphere::intersect(const Ray& ray, Hit& hit) const {
    // Treat as ellipsoid: p_world = C + S * p_local, with unit sphere in local.
    const RayVec3 ro_w = ray.getPosition();
    const RayVec3 rd_w = ray.getDirection();

    // Transform to local by inverse scale around center
    const double invSx = (m_scale.x != 0.0) ? (1.0 / m_scale.x) : 0.0;
    const double invSy = (m_scale.y != 0.0) ? (1.0 / m_scale.y) : 0.0;
    const double invSz = (m_scale.z != 0.0) ? (1.0 / m_scale.z) : 0.0;

    const double ro_local[3] = { (ro_w.x - m_location.x) * invSx,
                                 (ro_w.y - m_location.y) * invSy,
                                 (ro_w.z - m_location.z) * invSz };
    const double rd_local[3] = { rd_w.x * invSx, rd_w.y * invSy, rd_w.z * invSz };

    // Intersect with sphere of radius r_local. If radius provided, scale unit sphere by r first.
    const double r = 1.0; // unit sphere in local metrics; world scale shapes it
    const double A = rd_local[0]*rd_local[0] + rd_local[1]*rd_local[1] + rd_local[2]*rd_local[2];
    const double B = 2.0 * (ro_local[0]*rd_local[0] + ro_local[1]*rd_local[1] + ro_local[2]*rd_local[2]);
    const double C = ro_local[0]*ro_local[0] + ro_local[1]*ro_local[1] + ro_local[2]*ro_local[2] - r*r;
    const double eps = 1e-12;
    if (A < eps) return false;
    const double disc = B*B - 4.0*A*C;
    if (disc < 0.0) return false;
    const double sqrt_disc = std::sqrt(disc);
    double t0 = (-B - sqrt_disc) / (2.0*A);
    double t1 = (-B + sqrt_disc) / (2.0*A);
    if (t0 > t1) std::swap(t0, t1);
    double t_local = (t0 >= 0.0) ? t0 : t1;
    if (t_local < 0.0) return false;

    // Local intersection point
    const double p_local[3] = { ro_local[0] + t_local*rd_local[0],
                                ro_local[1] + t_local*rd_local[1],
                                ro_local[2] + t_local*rd_local[2] };

    // Transform point back to world: p_w = C + S * p_local
    const double p_world[3] = { m_location.x + p_local[0] * m_scale.x,
                                m_location.y + p_local[1] * m_scale.y,
                                m_location.z + p_local[2] * m_scale.z };

    // Normal: n_world ∝ S^{-1} * n_local, where n_local = p_local / r
    double n_local[3] = { p_local[0] / r, p_local[1] / r, p_local[2] / r };
    double n_world_raw[3] = { n_local[0] * invSx, n_local[1] * invSy, n_local[2] * invSz };
    const double nlen = std::sqrt(n_world_raw[0]*n_world_raw[0] + n_world_raw[1]*n_world_raw[1] + n_world_raw[2]*n_world_raw[2]);
    double n_unit[3] = { n_world_raw[0]/nlen, n_world_raw[1]/nlen, n_world_raw[2]/nlen };

    // Reflection in world space
    const double d_world[3] = { rd_w.x, rd_w.y, rd_w.z };
    const double d_dot_n = d_world[0]*n_unit[0] + d_world[1]*n_unit[1] + d_world[2]*n_unit[2];
    double r_world_tmp[3] = { d_world[0] - 2.0 * d_dot_n * n_unit[0],
                              d_world[1] - 2.0 * d_dot_n * n_unit[1],
                              d_world[2] - 2.0 * d_dot_n * n_unit[2] };
    double rlen = std::sqrt(r_world_tmp[0]*r_world_tmp[0] + r_world_tmp[1]*r_world_tmp[1] + r_world_tmp[2]*r_world_tmp[2]);
    if (rlen < 1e-15) rlen = 1.0;
    double r_world[3] = { r_world_tmp[0]/rlen, r_world_tmp[1]/rlen, r_world_tmp[2]/rlen };

    HitVec3 ip{ p_world[0], p_world[1], p_world[2] };
    HitVec3 nrm{ n_unit[0], n_unit[1], n_unit[2] };
    HitVec3 refl{ r_world[0], r_world[1], r_world[2] };
    hit.setIntersectionPoint(ip);
    hit.setSurfaceNormal(nrm);
    hit.setReflectedDirection(refl);

    // Distance along original world ray: project vector from origin to ip onto direction
    const double dx = p_world[0] - ro_w.x;
    const double dy = p_world[1] - ro_w.y;
    const double dz = p_world[2] - ro_w.z;
    const double d_len = std::sqrt(d_world[0]*d_world[0] + d_world[1]*d_world[1] + d_world[2]*d_world[2]);
    const double dist_world = (d_len > 1e-15) ? (dx*d_world[0] + dy*d_world[1] + dz*d_world[2]) / d_len
                                              : std::sqrt(dx*dx + dy*dy + dz*dz);
    hit.setDistanceAlongRay(dist_world);
    return true;
}

void Sphere::compute_aabb(Float3& outMin, Float3& outMax) const {
    // Transform local cube [-1,1]^3 by non-uniform scale and translation
    outMin = { +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity() };
    outMax = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
    auto transform_point = [&](double lx, double ly, double lz, double out[3]) {
        out[0] = m_location.x + lx * m_scale.x;
        out[1] = m_location.y + ly * m_scale.y;
        out[2] = m_location.z + lz * m_scale.z;
    };
    for (int dx = -1; dx <= 1; dx += 2) {
        for (int dy = -1; dy <= 1; dy += 2) {
            for (int dz = -1; dz <= 1; dz += 2) {
                double pw[3];
                transform_point(static_cast<double>(dx), static_cast<double>(dy), static_cast<double>(dz), pw);
                aabb_expand(outMin, outMax, pw);
            }
        }
    }
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


bool Plane::intersect(const Ray& ray, Hit& hit) const {
    if (m_corners.size() < 4) return false;

    // Build local basis from corners: P(u,v,w) = P0 + u*e1 + v*e2 + w*n, with plane at w=0
    const Float3& P0 = m_corners[0];
    const Float3& P1 = m_corners[1];
    const Float3& P3 = m_corners[3];

    const double e1[3] = { P1.x - P0.x, P1.y - P0.y, P1.z - P0.z };
    const double e2[3] = { P3.x - P0.x, P3.y - P0.y, P3.z - P0.z };

    auto cross3 = [](const double a[3], const double b[3], double out[3]) {
        out[0] = a[1]*b[2] - a[2]*b[1];
        out[1] = a[2]*b[0] - a[0]*b[2];
        out[2] = a[0]*b[1] - a[1]*b[0];
    };
    auto dot3 = [](const double a[3], const double b[3]) -> double {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };

    double n_col[3];
    cross3(e1, e2, n_col);
    const double n_len = std::sqrt(n_col[0]*n_col[0] + n_col[1]*n_col[1] + n_col[2]*n_col[2]);
    if (n_len < 1e-15) return false; // degenerate
    const double n_unit[3] = { n_col[0]/n_len, n_col[1]/n_len, n_col[2]/n_len };

    // Matrix M = [e1 e2 n] and its inverse to transform world -> local
    const double M[3][3] = {
        { e1[0], e2[0], n_unit[0] },
        { e1[1], e2[1], n_unit[1] },
        { e1[2], e2[2], n_unit[2] }
    };

    auto det3 = [&](const double A[3][3]) -> double {
        return A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
             - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
             + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    };
    const double detM = det3(M);
    if (std::abs(detM) < 1e-15) return false;
    const double invDet = 1.0 / detM;

    double Minv[3][3];
    // Adjugate (cofactor transpose) divided by det
    Minv[0][0] =  (M[1][1]*M[2][2] - M[1][2]*M[2][1]) * invDet;
    Minv[0][1] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) * invDet;
    Minv[0][2] =  (M[0][1]*M[1][2] - M[0][2]*M[1][1]) * invDet;
    Minv[1][0] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) * invDet;
    Minv[1][1] =  (M[0][0]*M[2][2] - M[0][2]*M[2][0]) * invDet;
    Minv[1][2] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) * invDet;
    Minv[2][0] =  (M[1][0]*M[2][1] - M[1][1]*M[2][0]) * invDet;
    Minv[2][1] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) * invDet;
    Minv[2][2] =  (M[0][0]*M[1][1] - M[0][1]*M[1][0]) * invDet;

    auto mul3 = [](const double A[3][3], const double v[3], double out[3]) {
        out[0] = A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
        out[1] = A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
        out[2] = A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2];
    };

    // World ray
    const RayVec3 ro_w = ray.getPosition();
    const RayVec3 rd_w = ray.getDirection();

    // Transform to local coordinates of the plane
    const double ro_rel[3] = { ro_w.x - P0.x, ro_w.y - P0.y, ro_w.z - P0.z };
    double ro_l[3]; mul3(Minv, ro_rel, ro_l);
    const double rd_in[3] = { rd_w.x, rd_w.y, rd_w.z };
    double rd_l[3]; mul3(Minv, rd_in, rd_l);

    // Intersect local plane w-axis is z, plane is z=0
    const double eps = 1e-12;
    if (std::abs(rd_l[2]) < eps) {
        if (std::abs(ro_l[2]) >= eps) return false; // parallel and not on plane
        return false; // treat as no unique intersection
    }

    const double t_hit = -ro_l[2] / rd_l[2];
    if (t_hit < 0.0) return false;

    const double u = ro_l[0] + t_hit * rd_l[0];
    const double v = ro_l[1] + t_hit * rd_l[1];

    // Inside quad if 0<=u<=1 and 0<=v<=1 (with small epsilon tolerance)
    if (u < -eps || u > 1.0 + eps || v < -eps || v > 1.0 + eps) return false;

    // Local intersection point
    const double p_l[3] = { u, v, 0.0 };

    // Transform back to world: p_w = P0 + M * p_l
    double p_w_rel[3]; mul3(M, p_l, p_w_rel);
    double p_world[3] = { P0.x + p_w_rel[0], P0.y + p_w_rel[1], P0.z + p_w_rel[2] };

    // World-space normal (already unit length)
    const double n_world[3] = { n_unit[0], n_unit[1], n_unit[2] };

    // Reflected direction in world space
    const double d_world[3] = { rd_w.x, rd_w.y, rd_w.z };
    const double d_dot_n = dot3(d_world, n_world);
    double r_world_tmp[3] = { d_world[0] - 2.0 * d_dot_n * n_world[0],
                              d_world[1] - 2.0 * d_dot_n * n_world[1],
                              d_world[2] - 2.0 * d_dot_n * n_world[2] };
    double r_len = std::sqrt(r_world_tmp[0]*r_world_tmp[0] + r_world_tmp[1]*r_world_tmp[1] + r_world_tmp[2]*r_world_tmp[2]);
    if (r_len < 1e-15) r_len = 1.0;
    double r_world[3] = { r_world_tmp[0]/r_len, r_world_tmp[1]/r_len, r_world_tmp[2]/r_len };

    // Populate Hit
    HitVec3 ip{ p_world[0], p_world[1], p_world[2] };
    HitVec3 nrm{ n_world[0], n_world[1], n_world[2] };
    HitVec3 refl{ r_world[0], r_world[1], r_world[2] };

    hit.setIntersectionPoint(ip);
    hit.setSurfaceNormal(nrm);
    hit.setReflectedDirection(refl);

    const double dx = p_world[0] - ro_w.x;
    const double dy = p_world[1] - ro_w.y;
    const double dz = p_world[2] - ro_w.z;
    const double d_len = std::sqrt(d_world[0]*d_world[0] + d_world[1]*d_world[1] + d_world[2]*d_world[2]);
    const double dist_world = (d_len > 1e-15) ? (dx*d_world[0] + dy*d_world[1] + dz*d_world[2]) / d_len
                                              : std::sqrt(dx*dx + dy*dy + dz*dz);
    hit.setDistanceAlongRay(dist_world);
    return true;
}

void Plane::compute_aabb(Float3& outMin, Float3& outMax) const {
    if (m_corners.empty()) { outMin = {0,0,0}; outMax = {0,0,0}; return; }
    outMin = { +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity() };
    outMax = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
    for (const auto& c : m_corners) {
        double p[3] = { c.x, c.y, c.z };
        aabb_expand(outMin, outMax, p);
    }
}

// ---------- BoundingBox ----------

BoundingBox::BoundingBox() : m_min{0,0,0}, m_max{0,0,0}, m_left(nullptr), m_right(nullptr) {}
BoundingBox::BoundingBox(const Float3& bmin, const Float3& bmax) : m_min(bmin), m_max(bmax), m_left(nullptr), m_right(nullptr) {}

void BoundingBox::write_to_console(std::ostream& out) const {
    out << "BoundingBox(min=[" << m_min.x << ", " << m_min.y << ", " << m_min.z
        << "], max=[" << m_max.x << ", " << m_max.y << ", " << m_max.z << "])\n";
}

bool BoundingBox::intersect(const Ray& ray, Hit& hit) const {
    const RayVec3 ro = ray.getPosition();
    const RayVec3 rd = ray.getDirection();
    const double eps = 1e-12;
    double tmin = -std::numeric_limits<double>::infinity();
    double tmax =  std::numeric_limits<double>::infinity();
    auto slab = [&](double o, double d, double lo, double hi) {
        if (std::abs(d) < eps) {
            if (o < lo || o > hi) return false;
            return true;
        }
        double t1 = (lo - o) / d;
        double t2 = (hi - o) / d;
        if (t1 > t2) std::swap(t1, t2);
        if (t1 > tmin) tmin = t1;
        if (t2 < tmax) tmax = t2;
        return tmin <= tmax;
    };
    if (!slab(ro.x, rd.x, m_min.x, m_max.x)) return false;
    if (!slab(ro.y, rd.y, m_min.y, m_max.y)) return false;
    if (!slab(ro.z, rd.z, m_min.z, m_max.z)) return false;
    if (tmax < 0.0) return false;
    const double t_hit = (tmin >= 0.0) ? tmin : tmax;
    hit.setDistanceAlongRay(t_hit);
    return true;
}

void BoundingBox::compute_aabb(Float3& outMin, Float3& outMax) const { outMin = m_min; outMax = m_max; }

void BoundingBox::setChildren(const Mesh* left, const Mesh* right) { m_left = left; m_right = right; }
const Mesh* BoundingBox::getLeft() const { return m_left; }
const Mesh* BoundingBox::getRight() const { return m_right; }
const Float3& BoundingBox::getMin() const { return m_min; }
const Float3& BoundingBox::getMax() const { return m_max; }
void BoundingBox::setBounds(const Float3& bmin, const Float3& bmax) { m_min = bmin; m_max = bmax; }


