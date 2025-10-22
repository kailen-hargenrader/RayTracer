#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>
#include <iosfwd>

// Forward declarations to avoid including full headers here
class Hit;
class Ray;
struct Float3; // forward declare for Mesh::compute_aabb

/** Base class for all mesh types. */
class Mesh {
public:
    virtual ~Mesh() = default;

    /** Write a human-readable representation to the provided stream. */
    virtual void write_to_console(std::ostream& out) const = 0;

    /** Intersect this mesh with a ray; populate Hit on success. */
    virtual bool intersect(const Ray& ray, Hit& hit) const { return false; }

    /** Compute an axis-aligned bounding box in world coordinates for this mesh. */
    virtual void compute_aabb(Float3& outMin, Float3& outMax) const = 0;

    /** Return this as a BoundingBox if applicable, else nullptr (avoids RTTI cost). */
    virtual const class BoundingBox* asBoundingBox() const { return nullptr; }
};

/** Simple 3D point/vector container used by mesh shapes. */
struct Float3 {
    double x;
    double y;
    double z;
};

/** Euler angles (roll, pitch, yaw in radians). */
struct EulerAngles {
    double roll;
    double pitch;
    double yaw;
};

/** Unit cube mesh with transform (translation, rotation) and per-axis scale. */
class Cube : public Mesh {
public:
    Cube();
    Cube(const Float3& translation, const EulerAngles& rotation, const Float3& scale);

    static std::vector<Cube> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;

    bool intersect(const Ray& ray, Hit& hit) const override;

    void compute_aabb(Float3& outMin, Float3& outMax) const override;

private:
    Float3 m_translation;
    EulerAngles m_rotation;
    Float3 m_scale;
};

/** Cylinder mesh defined by transform and per-axis scale only.
 *  Canonical local cylinder: x^2 + y^2 = 1, z in [-1, 1]
 *  Forward transform: p_world = T + R * (S * p_local)
 */
class Cylinder : public Mesh {
public:
    Cylinder();
    Cylinder(const Float3& translation, const EulerAngles& rotation, const Float3& scale);

    static std::vector<Cylinder> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;
    bool intersect(const Ray& ray, Hit& hit) const override;

    void compute_aabb(Float3& outMin, Float3& outMax) const override;

private:
    Float3 m_translation;
    EulerAngles m_rotation;
    Float3 m_scale; // per-axis
};

/** Sphere mesh defined by center location and per-axis scale only.
 *  Canonical local sphere: x^2 + y^2 + z^2 = 1; world ellipsoid via scale.
 *  Forward transform: p_world = C + S * p_local
 */
class Sphere : public Mesh {
public:
    Sphere();
    Sphere(const Float3& location, const Float3& scale);

    static std::vector<Sphere> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;
    bool intersect(const Ray& ray, Hit& hit) const override;

    void compute_aabb(Float3& outMin, Float3& outMax) const override;

private:
    Float3 m_location;
    Float3 m_scale; // per-axis
};

/** Plane mesh defined by four corner points. */
class Plane : public Mesh {
public:
    Plane();
    explicit Plane(const std::vector<Float3>& corners);

    static std::vector<Plane> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;
    bool intersect(const Ray& ray, Hit& hit) const override;

    void compute_aabb(Float3& outMin, Float3& outMax) const override;

private:
    std::vector<Float3> m_corners; // expected size 4
};

/** Axis-aligned bounding box mesh used as internal nodes in BVH. */
class BoundingBox : public Mesh {
public:
    BoundingBox();
    BoundingBox(const Float3& bmin, const Float3& bmax);

    void write_to_console(std::ostream& out) const override;
    bool intersect(const Ray& ray, Hit& hit) const override;
    void compute_aabb(Float3& outMin, Float3& outMax) const override;

    const BoundingBox* asBoundingBox() const override { return this; }

    void setChildren(const Mesh* left, const Mesh* right);
    const Mesh* getLeft() const;
    const Mesh* getRight() const;

    const Float3& getMin() const;
    const Float3& getMax() const;
    void setBounds(const Float3& bmin, const Float3& bmax);

private:
    Float3 m_min;
    Float3 m_max;
    const Mesh* m_left;
    const Mesh* m_right;
};

#endif // MESH_H


