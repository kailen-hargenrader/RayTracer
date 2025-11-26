#include "raytracer/core/Scene.h"
#include "raytracer/utils/Json.h"
#include "raytracer/geometry/Sphere.h"
#include "raytracer/geometry/Cube.h"
#include "raytracer/geometry/Cylinder.h"
#include "raytracer/geometry/Plane.h"
#include "raytracer/accel/AABB.h"
#include "raytracer/math/Mat4.h"
#include <array>
#include <unordered_map>
#include <filesystem>

namespace rt {

static Material parseMaterial(const Json::Value* mat) {
    Material m;
    if (!mat || !mat->isObject()) return m;
    const auto& o = mat->asObject();
    auto itBase = o.find("baseColor");
    if (itBase != o.end() && itBase->second.isObject()) {
        const auto& bc = itBase->second.asObject();
        auto itType = bc.find("type");
        if (itType != bc.end() && itType->second.isString()) {
            if (itType->second.asString() == "RGB") {
                auto itVal = bc.find("value");
                if (itVal != bc.end() && itVal->second.isArray()) {
                    const auto& arr = itVal->second.asArray();
                    if (arr.size() >= 3 && arr[0].isNumber() && arr[1].isNumber() && arr[2].isNumber()) {
                        m.baseColor = Vec3(static_cast<float>(arr[0].asNumber()), static_cast<float>(arr[1].asNumber()), static_cast<float>(arr[2].asNumber()));
                    }
                }
            } else if (itType->second.asString() == "ImageTexture") {
                auto itPath = bc.find("path");
                if (itPath != bc.end() && itPath->second.isString()) m.baseColorTexturePath = itPath->second.asString();
            }
        }
    }
    auto f = [&](const char* k, float& dst){
        auto it = o.find(k);
        if (it != o.end() && it->second.isNumber()) dst = static_cast<float>(it->second.asNumber());
    };
    f("alpha", m.alpha);
    f("metallic", m.metallic);
    f("roughness", m.roughness);
    f("ior", m.ior);
    return m;
}

std::shared_ptr<Scene> Scene::loadFromJsonFile(const std::string& path) {
    Json::Value root;
    std::string err;
    if (!Json::parseFile(path, root, &err)) {
        return nullptr;
    }

    auto scene = std::make_shared<Scene>();
    if (!root.isObject()) return scene;
    const auto& ro = root.asObject();

    // Camera
    if (const Json::Value* camSec = root.get("Camera")) {
        if (camSec->isObject()) {
            const auto& co = camSec->asObject();
            auto it = co.find("Perspective");
            if (it != co.end() && it->second.isArray()) {
                const auto& arr = it->second.asArray();
                if (!arr.empty() && arr[0].isObject()) {
                    const auto& c0 = arr[0].asObject();
                    auto getVec3Key = [](const std::unordered_map<std::string, Json::Value>& obj, const char* key)->Vec3{
                        auto it = obj.find(key);
                        if (it == obj.end() || !it->second.isArray()) return {0,0,1};
                        const auto& a = it->second.asArray();
                        if (a.size()>=3) {
                            return { static_cast<float>(a[0].asNumber()), static_cast<float>(a[1].asNumber()), static_cast<float>(a[2].asNumber()) };
                        }
                        return {0,0,1};
                    };
                    Vec3 loc = getVec3Key(c0, "location");
                    Vec3 dir(0.0f, 0.0f, -1.0f); // default camera forward is -Z in local space
                    // Prefer explicit direction if provided
                    auto itDir = c0.find("direction");
                    if (itDir != c0.end() && itDir->second.isArray() && itDir->second.asArray().size() >= 3) {
                        const auto& a = itDir->second.asArray();
                        dir = { static_cast<float>(a[0].asNumber()), static_cast<float>(a[1].asNumber()), static_cast<float>(a[2].asNumber()) };
                    } else {
                        // Fallback to rotationEuler if present
                        auto itRot = c0.find("rotationEuler");
                        bool hasRot = (itRot != c0.end()) && itRot->second.isArray() && (itRot->second.asArray().size() >= 3);
                        if (hasRot) {
                            const auto& ra = itRot->second.asArray();
                            float rx = static_cast<float>(ra[0].asNumber());
                            float ry = static_cast<float>(ra[1].asNumber());
                            float rz = static_cast<float>(ra[2].asNumber());
                            Mat4 R = Mat4::rotationXYZ(rx, ry, rz);
                            // Assume camera looks along -Z in its local space
                            dir = R.transformDirection({0.0f, 0.0f, -1.0f});
                        }
                    }
                    dir = dir.normalized();
                    float focal = 35.0f, sensorW = 36.0f, sensorH = 24.0f;
                    auto fnum = [&](const char* k, float& dst){
                        auto itv = c0.find(k);
                        if (itv != c0.end() && itv->second.isNumber()) dst = static_cast<float>(itv->second.asNumber());
                    };
                    fnum("focalLengthMM", focal);
                    fnum("sensorWidthMM", sensorW);
                    fnum("sensorHeightMM", sensorH);
                    int resW = 640, resH = 480;
                    auto itFR = c0.find("filmResolution");
                    if (itFR != c0.end()) {
                        const Json::Value* fr = &itFR->second;
                        if (fr->isArray() && fr->asArray().size()>=2) {
                            resW = static_cast<int>(fr->asArray()[0].asNumber());
                            resH = static_cast<int>(fr->asArray()[1].asNumber());
                        }
                    }
                    scene->camera.set(loc, dir, focal, sensorW, sensorH, resW, resH);
                }
            }
        }
    }

    // Lights
    if (const Json::Value* lightSec = root.get("Light")) {
        if (lightSec->isObject()) {
            const auto& lo = lightSec->asObject();
            auto it = lo.find("Point");
            if (it != lo.end() && it->second.isArray()) {
                for (const auto& v : it->second.asArray()) {
                    if (!v.isObject()) continue;
                    const auto& po = v.asObject();
                    PointLight L;
                    if (const Json::Value* loc = v.get("location")) {
                        if (loc->isArray() && loc->asArray().size()>=3) {
                            L.position = { static_cast<float>(loc->asArray()[0].asNumber()),
                                           static_cast<float>(loc->asArray()[1].asNumber()),
                                           static_cast<float>(loc->asArray()[2].asNumber()) };
                        }
                    }
                    if (const Json::Value* ri = v.get("radiantIntensity")) {
                        L.intensity = static_cast<float>(ri->asNumber());
                    } else {
                        L.intensity = 10.0f;
                    }
                    scene->lights.push_back(L);
                }
            }
        }
    }

    // Meshes
    if (const Json::Value* meshSec = root.get("Mesh")) {
        if (meshSec->isObject()) {
            auto readTRS = [](const Json::Value& o, Vec3& t, Vec3& r, Vec3& s){
                auto arrToVec3 = [](const Json::Value* v)->Vec3{
                    if (!v || !v->isArray() || v->asArray().size()<3) return {0,0,0};
                    return { static_cast<float>(v->asArray()[0].asNumber()),
                             static_cast<float>(v->asArray()[1].asNumber()),
                             static_cast<float>(v->asArray()[2].asNumber()) };
                };
                t = arrToVec3(o.get("translation"));
                r = arrToVec3(o.get("rotationEuler"));
                s = arrToVec3(o.get("scale"));
            };

            auto addObj = [&](const std::string& kind, const Json::Value& v){
                if (!v.isArray()) return;
                for (const auto& it : v.asArray()) {
                    if (!it.isObject()) continue;
                    const auto& oo = it.asObject();
                    Material mat = parseMaterial(oo.count("material") ? &oo.at("material") : nullptr);

                    if (kind == "Plane") {
                        // corners world-space
                        const Json::Value* corners = it.get("corners");
                        if (!corners || !corners->isArray() || corners->asArray().size()<4) continue;
                        std::array<Vec3,4> cs;
                        for (int i = 0; i < 4; ++i) {
                            const auto& c = corners->asArray()[i];
                            cs[i] = { static_cast<float>(c.asArray()[0].asNumber()),
                                      static_cast<float>(c.asArray()[1].asNumber()),
                                      static_cast<float>(c.asArray()[2].asNumber()) };
                        }
                        auto pl = std::make_shared<Plane>();
                        pl->setCorners(cs);
                        SceneObject so;
                        so.mesh = pl;
                        so.material = mat;
                        // Resolve texture absolute path if present
                        if (!so.material.baseColorTexturePath.empty()) {
                            std::filesystem::path jsonPath(path);
                            auto base = jsonPath.parent_path().parent_path(); // s2899221/
                            std::filesystem::path texAbs = base / std::filesystem::path(so.material.baseColorTexturePath);
                            so.material.baseColorTexturePath = texAbs.string();
                        }
                        // Bounds from corners
                        AABB b;
                        for (const auto& p : cs) b.expand(p);
                        so.bounds = b;
                        scene->objects.push_back(std::move(so));
                    } else {
                        Vec3 t, r, s;
                        readTRS(it, t, r, s);
                        // Blender primitives (Cube/Sphere/Cylinder) default to size 2 units (side=2, radius=1, height=2).
                        // Our canonical primitives use unit size with half-extents/radius = 0.5.
                        // Multiply incoming scales by 2 to match Blender's apparent size.
                        s = s * 2.0f;
                        // Optional quaternion rotation (world-space), preferred over Euler if present
                        float qw = 0.0f, qx = 0.0f, qy = 0.0f, qz = 0.0f;
                        bool hasQuat = false;
                        if (const Json::Value* rq = it.get("rotationQuat")) {
                            if (rq->isArray() && rq->asArray().size() >= 4) {
                                const auto& a = rq->asArray();
                                if (a[0].isNumber() && a[1].isNumber() && a[2].isNumber() && a[3].isNumber()) {
                                    qw = static_cast<float>(a[0].asNumber());
                                    qx = static_cast<float>(a[1].asNumber());
                                    qy = static_cast<float>(a[2].asNumber());
                                    qz = static_cast<float>(a[3].asNumber());
                                    hasQuat = true;
                                }
                            }
                        }
                        Mat4 M;
                        if (hasQuat) {
                            M = Mat4::translation(t) * Mat4::rotationQuaternion(qw, qx, qy, qz) * Mat4::scaling(s);
                        } else {
                            M = Mat4::TRS(t, r, s);
                        }
                        MeshPtr mp;
                        if (kind == "Sphere") {
                            auto sp = std::make_shared<Sphere>();
                            sp->setTransform(M);
                            mp = sp;
                        } else if (kind == "Cube") {
                            auto cu = std::make_shared<Cube>();
                            cu->setTransform(M);
                            mp = cu;
                        } else if (kind == "Cylinder") {
                            auto cy = std::make_shared<Cylinder>();
                            cy->setTransform(M);
                            mp = cy;
                        } else {
                            continue;
                        }
                        SceneObject so;
                        so.mesh = mp;
                        so.material = mat;
                        // Transform canonical bounds using TRS
                        // We cannot access m_objectToWorld here (protected). Reconstruct it with TRS:
                        so.bounds = transformAABB(AABB({-0.5f,-0.5f,-0.5f},{0.5f,0.5f,0.5f}), M);
                        // Resolve texture absolute path if present
                        if (!so.material.baseColorTexturePath.empty()) {
                            std::filesystem::path jsonPath(path);
                            auto base = jsonPath.parent_path().parent_path(); // s2899221/
                            std::filesystem::path texAbs = base / std::filesystem::path(so.material.baseColorTexturePath);
                            so.material.baseColorTexturePath = texAbs.string();
                        }
                        scene->objects.push_back(std::move(so));
                    }
                }
            };

            const auto& mo = meshSec->asObject();
            for (const auto& kv : mo) {
                addObj(kv.first, kv.second);
            }
        }
    }

    return scene;
}

} // namespace rt


