#include "raytracer/math/Mat4.h"
#include <cstring>
#include <algorithm>

namespace rt {

Mat4::Mat4() {
    std::memset(m, 0, sizeof(m));
}

Mat4 Mat4::identity() {
    Mat4 I;
    I.m[0][0] = 1.0f; I.m[1][1] = 1.0f; I.m[2][2] = 1.0f; I.m[3][3] = 1.0f;
    return I;
}

Mat4 Mat4::translation(const Vec3& t) {
    Mat4 M = identity();
    M.m[0][3] = t.x;
    M.m[1][3] = t.y;
    M.m[2][3] = t.z;
    return M;
}

Mat4 Mat4::scaling(const Vec3& s) {
    Mat4 M = identity();
    M.m[0][0] = s.x;
    M.m[1][1] = s.y;
    M.m[2][2] = s.z;
    return M;
}

Mat4 Mat4::rotationXYZ(float rx, float ry, float rz) {
    // Rotation order: X then Y then Z (Rz * Ry * Rx when applied to column vectors)
    float cx = std::cos(rx), sx = std::sin(rx);
    float cy = std::cos(ry), sy = std::sin(ry);
    float cz = std::cos(rz), sz = std::sin(rz);

    Mat4 Rx = identity();
    Rx.m[1][1] = cx; Rx.m[1][2] = -sx;
    Rx.m[2][1] = sx; Rx.m[2][2] = cx;

    Mat4 Ry = identity();
    Ry.m[0][0] = cy; Ry.m[0][2] = sy;
    Ry.m[2][0] = -sy; Ry.m[2][2] = cy;

    Mat4 Rz = identity();
    Rz.m[0][0] = cz; Rz.m[0][1] = -sz;
    Rz.m[1][0] = sz; Rz.m[1][1] = cz;

    // Combined: R = Rz * Ry * Rx
    return multiply(multiply(Rz, Ry), Rx);
}

Mat4 Mat4::rotationQuaternion(float w, float x, float y, float z) {
    // Normalize quaternion to avoid scaling the rotation
    float n = std::sqrt(w*w + x*x + y*y + z*z);
    if (n > 1e-8f) {
        float inv = 1.0f / n;
        w *= inv; x *= inv; y *= inv; z *= inv;
    } else {
        return identity();
    }

    // Column-vector convention, row-major storage
    // Rotation matrix derived from quaternion (w, x, y, z)
    float xx = x * x, yy = y * y, zz = z * z;
    float xy = x * y, xz = x * z, yz = y * z;
    float wx = w * x, wy = w * y, wz = w * z;

    Mat4 R = identity();
    R.m[0][0] = 1.0f - 2.0f * (yy + zz);
    R.m[0][1] = 2.0f * (xy - wz);
    R.m[0][2] = 2.0f * (xz + wy);

    R.m[1][0] = 2.0f * (xy + wz);
    R.m[1][1] = 1.0f - 2.0f * (xx + zz);
    R.m[1][2] = 2.0f * (yz - wx);

    R.m[2][0] = 2.0f * (xz - wy);
    R.m[2][1] = 2.0f * (yz + wx);
    R.m[2][2] = 1.0f - 2.0f * (xx + yy);
    return R;
}

Mat4 Mat4::TRS(const Vec3& t, const Vec3& rEuler, const Vec3& s) {
    return translation(t) * rotationXYZ(rEuler.x, rEuler.y, rEuler.z) * scaling(s);
}

Mat4 Mat4::multiply(const Mat4& a, const Mat4& b) {
    Mat4 r;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            r.m[i][j] = a.m[i][0] * b.m[0][j] +
                        a.m[i][1] * b.m[1][j] +
                        a.m[i][2] * b.m[2][j] +
                        a.m[i][3] * b.m[3][j];
        }
    }
    return r;
}

Mat4 Mat4::inverseAffine() const {
    // Assumes last row = [0 0 0 1]
    // Inverse of [ R t ; 0 1 ] is [ R^{-1}  -R^{-1} t ; 0 1 ]
    // Where R is upper-left 3x3 (including non-uniform scaling and rotation)
    // We'll compute inverse of 3x3 using adjugate/determinant.
    const float a00 = m[0][0], a01 = m[0][1], a02 = m[0][2];
    const float a10 = m[1][0], a11 = m[1][1], a12 = m[1][2];
    const float a20 = m[2][0], a21 = m[2][1], a22 = m[2][2];

    float det = a00*(a11*a22 - a12*a21) - a01*(a10*a22 - a12*a20) + a02*(a10*a21 - a11*a20);
    if (std::abs(det) < 1e-8f) {
        // Fallback to identity to avoid NaNs
        return identity();
    }
    float invDet = 1.0f / det;

    Mat4 inv = identity();
    // Inverse 3x3 (cofactor transpose) * 1/det
    inv.m[0][0] =  (a11*a22 - a12*a21) * invDet;
    inv.m[0][1] = -(a01*a22 - a02*a21) * invDet;
    inv.m[0][2] =  (a01*a12 - a02*a11) * invDet;

    inv.m[1][0] = -(a10*a22 - a12*a20) * invDet;
    inv.m[1][1] =  (a00*a22 - a02*a20) * invDet;
    inv.m[1][2] = -(a00*a12 - a02*a10) * invDet;

    inv.m[2][0] =  (a10*a21 - a11*a20) * invDet;
    inv.m[2][1] = -(a00*a21 - a01*a20) * invDet;
    inv.m[2][2] =  (a00*a11 - a01*a10) * invDet;

    // Translation
    const float tx = m[0][3], ty = m[1][3], tz = m[2][3];
    // -R^{-1} t
    inv.m[0][3] = -(inv.m[0][0]*tx + inv.m[0][1]*ty + inv.m[0][2]*tz);
    inv.m[1][3] = -(inv.m[1][0]*tx + inv.m[1][1]*ty + inv.m[1][2]*tz);
    inv.m[2][3] = -(inv.m[2][0]*tx + inv.m[2][1]*ty + inv.m[2][2]*tz);
    return inv;
}

Vec3 Mat4::transformPoint(const Vec3& p) const {
    float x = m[0][0]*p.x + m[0][1]*p.y + m[0][2]*p.z + m[0][3];
    float y = m[1][0]*p.x + m[1][1]*p.y + m[1][2]*p.z + m[1][3];
    float z = m[2][0]*p.x + m[2][1]*p.y + m[2][2]*p.z + m[2][3];
    float w = m[3][0]*p.x + m[3][1]*p.y + m[3][2]*p.z + m[3][3];
    if (std::abs(w) > 1e-8f) {
        float iw = 1.0f / w;
        return {x * iw, y * iw, z * iw};
    }
    return {x, y, z};
}

Vec3 Mat4::transformVector(const Vec3& v) const {
    // Applies translation
    float x = m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z + m[0][3];
    float y = m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z + m[1][3];
    float z = m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z + m[2][3];
    return {x, y, z};
}

Vec3 Mat4::transformDirection(const Vec3& d) const {
    // Ignores translation
    float x = m[0][0]*d.x + m[0][1]*d.y + m[0][2]*d.z;
    float y = m[1][0]*d.x + m[1][1]*d.y + m[1][2]*d.z;
    float z = m[2][0]*d.x + m[2][1]*d.y + m[2][2]*d.z;
    return {x, y, z};
}

Vec3 Mat4::transformNormal(const Vec3& n) const {
    // Using inverse transpose of upper-left 3x3 from this matrix
    // For our usage, we want to transform normals by (inverse(M)).T
    // This method expects that this Mat4 is the inverse-affine of the object-to-world matrix,
    // so that we can apply its transpose to normals.
    // To keep it simple, we'll compute it directly here using the upper-left 3x3.
    float a00 = m[0][0], a01 = m[1][0], a02 = m[2][0];
    float a10 = m[0][1], a11 = m[1][1], a12 = m[2][1];
    float a20 = m[0][2], a21 = m[1][2], a22 = m[2][2];
    Vec3 r(
        a00 * n.x + a01 * n.y + a02 * n.z,
        a10 * n.x + a11 * n.y + a12 * n.z,
        a20 * n.x + a21 * n.y + a22 * n.z
    );
    return r;
}

} // namespace rt





