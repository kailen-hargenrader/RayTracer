#pragma once

#include "raytracer/math/Vec3.h"
#include "raytracer/utils/Ray.h"
#include <cmath>

namespace rt {

class Camera {
public:
    // Set via Blender export data
    void set(const Vec3& position, const Vec3& direction, float focalLengthMM, float sensorWidthMM, float sensorHeightMM, int pixelWidth, int pixelHeight) {
        m_position = position;
        m_forward = direction.normalized();
        // Blender is Z-up; prefer Z as world up and fall back to Y if parallel
        Vec3 upWorld(0.0f, 0.0f, 1.0f);
        if (std::abs(Vec3::dot(upWorld, m_forward)) > 0.999f) upWorld = {0.0f, 1.0f, 0.0f};
        m_right = Vec3::cross(m_forward, upWorld).normalized();
        m_up = Vec3::cross(m_right, m_forward).normalized();

        m_pixelWidth = pixelWidth;
        m_pixelHeight = pixelHeight;
        // Compute FOV consistent with Blender's Sensor Fit = AUTO (approximation):
        // AUTO selects HORIZONTAL if image width >= image height, else VERTICAL.
        float imgAspect = (pixelHeight > 0) ? (static_cast<float>(pixelWidth) / static_cast<float>(pixelHeight)) : 1.0f;
        bool horizontalFit = pixelWidth >= pixelHeight;
        if (horizontalFit) {
            m_fovX = 2.0f * std::atan((0.5f * sensorWidthMM) / focalLengthMM);
            float halfTanX = std::tan(0.5f * m_fovX);
            m_fovY = 2.0f * std::atan(halfTanX / std::max(1e-6f, imgAspect));
        } else {
            m_fovY = 2.0f * std::atan((0.5f * sensorHeightMM) / focalLengthMM);
            float halfTanY = std::tan(0.5f * m_fovY);
            m_fovX = 2.0f * std::atan(halfTanY * imgAspect);
        }
    }

    int width() const { return m_pixelWidth; }
    int height() const { return m_pixelHeight; }

    // Sample in [0,1)^2 jitter space
    Ray generateRay(int px, int py, float jitterX = 0.5f, float jitterY = 0.5f) const {
        float ndcX = (static_cast<float>(px) + jitterX) / static_cast<float>(m_pixelWidth);
        float ndcY = (static_cast<float>(py) + jitterY) / static_cast<float>(m_pixelHeight);
        float screenX = (2.0f * ndcX - 1.0f) * std::tan(0.5f * m_fovX);
        float screenY = (1.0f - 2.0f * ndcY) * std::tan(0.5f * m_fovY);
        Vec3 dir = (m_forward + screenX * m_right + screenY * m_up).normalized();
        return Ray(m_position, dir);
    }

private:
    Vec3 m_position;
    Vec3 m_forward;
    Vec3 m_right;
    Vec3 m_up;
    float m_fovX {60.0f * 3.1415926535f / 180.0f};
    float m_fovY {40.0f * 3.1415926535f / 180.0f};
    int m_pixelWidth {640};
    int m_pixelHeight {480};
};

} // namespace rt


