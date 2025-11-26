#pragma once

#include <cmath>
#include <algorithm>
#include <cfloat>

namespace rt {

struct Vec3 {
    float x {0.0f};
    float y {0.0f};
    float z {0.0f};

    Vec3() = default;
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    float& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    float operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(float s) const { return {x / s, y / s, z / s}; }
    Vec3 operator*(const Vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }
    Vec3& operator/=(float s) { x /= s; y /= s; z /= s; return *this; }

    float length() const { return std::sqrt(x*x + y*y + z*z); }
    float lengthSquared() const { return x*x + y*y + z*z; }
    Vec3 normalized() const {
        float len = length();
        return (len > 0.0f) ? (*this / len) : *this;
    }

    static float dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
    static Vec3 cross(const Vec3& a, const Vec3& b) {
        return Vec3(
            a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x
        );
    }
    static Vec3 min(const Vec3& a, const Vec3& b) { return {std::min(a.x,b.x), std::min(a.y,b.y), std::min(a.z,b.z)}; }
    static Vec3 max(const Vec3& a, const Vec3& b) { return {std::max(a.x,b.x), std::max(a.y,b.y), std::max(a.z,b.z)}; }
    static Vec3 clamp01(const Vec3& v) {
        auto clamp = [](float t){ return std::max(0.0f, std::min(1.0f, t)); };
        return { clamp(v.x), clamp(v.y), clamp(v.z) };
    }
};

inline Vec3 operator*(float s, const Vec3& v) { return v * s; }

} // namespace rt


