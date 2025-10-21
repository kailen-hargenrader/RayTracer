#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>
#include <iosfwd>

// Forward declarations to avoid including full headers here
class Hit;
class Ray;

/** Base class for all mesh types. */
class Mesh {
public:
    virtual ~Mesh() = default;

    /** Write a human-readable representation to the provided stream. */
    virtual void write_to_console(std::ostream& out) const = 0;

    /** Intersect this mesh with a ray; populate Hit on success. */
    virtual bool intersect(const Ray& ray, Hit& hit) const { return false; }
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

/** Unit cube mesh with transform (translation, rotation) and uniform scale. */
class Cube : public Mesh {
public:
    Cube();
    Cube(const Float3& translation, const EulerAngles& rotation, double scale);

    static std::vector<Cube> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;

    bool intersect(const Ray& ray, Hit& hit) const override;

private:
    Float3 m_translation;
    EulerAngles m_rotation;
    double m_scale;
};

/** Cylinder mesh defined by transform, uniform scale, radius, and length. */
class Cylinder : public Mesh {
public:
    Cylinder();
    Cylinder(const Float3& translation, const EulerAngles& rotation, double scale, double radius, double length);

    static std::vector<Cylinder> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;

private:
    Float3 m_translation;
    EulerAngles m_rotation;
    double m_scale;
    double m_radius;
    double m_length;
};

/** Sphere mesh defined by center location and radius. */
class Sphere : public Mesh {
public:
    Sphere();
    Sphere(const Float3& location, double radius);

    static std::vector<Sphere> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;

private:
    Float3 m_location;
    double m_radius;
};

/** Plane mesh defined by four corner points. */
class Plane : public Mesh {
public:
    Plane();
    explicit Plane(const std::vector<Float3>& corners);

    static std::vector<Plane> read_from_json(const std::string& class_block);
    void write_to_console(std::ostream& out) const override;

private:
    std::vector<Float3> m_corners; // expected size 4
};

#endif // MESH_H


