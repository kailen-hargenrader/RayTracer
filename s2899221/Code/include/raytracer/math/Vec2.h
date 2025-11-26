#pragma once

#include <cmath>
#include <algorithm>

namespace rt {

struct Vec2 {
    float x {0.0f};
    float y {0.0f};

    Vec2() = default;
    Vec2(float x_, float y_) : x(x_), y(y_) {}

    float& operator[](int i) { return i == 0 ? x : y; }
    float operator[](int i) const { return i == 0 ? x : y; }

    Vec2 operator+(const Vec2& v) const { return {x + v.x, y + v.y}; }
    Vec2 operator-(const Vec2& v) const { return {x - v.x, y - v.y}; }
    Vec2 operator*(float s) const { return {x * s, y * s}; }
    Vec2 operator/(float s) const { return {x / s, y / s}; }
    Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    Vec2& operator*=(float s) { x *= s; y *= s; return *this; }
    Vec2& operator/=(float s) { x /= s; y /= s; return *this; }

    float length() const { return std::sqrt(x*x + y*y); }
    float lengthSquared() const { return x*x + y*y; }
    Vec2 normalized() const {
        float len = length();
        return (len > 0.0f) ? (*this / len) : *this;
    }
    static float dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
};

inline Vec2 operator*(float s, const Vec2& v) { return v * s; }

} // namespace rt


