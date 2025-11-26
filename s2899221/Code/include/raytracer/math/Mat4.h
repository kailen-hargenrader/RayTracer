#pragma once

#include <array>
#include <cmath>
#include "raytracer/math/Vec3.h"

namespace rt {

struct Mat4 {
    // Row-major 4x4
    float m[4][4];

    Mat4();
    static Mat4 identity();
    static Mat4 translation(const Vec3& t);
    static Mat4 scaling(const Vec3& s);
    static Mat4 rotationXYZ(float rx, float ry, float rz); // radians, XYZ order
    static Mat4 rotationQuaternion(float w, float x, float y, float z); // normalized quaternion
    static Mat4 TRS(const Vec3& t, const Vec3& rEuler, const Vec3& s);

    static Mat4 multiply(const Mat4& a, const Mat4& b);
    Mat4 operator*(const Mat4& b) const { return multiply(*this, b); }

    // For affine transforms (last row == [0,0,0,1])
    Mat4 inverseAffine() const;

    Vec3 transformPoint(const Vec3& p) const;
    Vec3 transformVector(const Vec3& v) const;   // includes translation
    Vec3 transformDirection(const Vec3& d) const; // ignores translation

    // 3x3 normal matrix derived from inverse transpose of upper-left
    Vec3 transformNormal(const Vec3& n) const;
};

} // namespace rt


