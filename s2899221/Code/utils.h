#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <optional>
#include <stdexcept>
#include <vector>

namespace util_json {

// Extracts the inner contents of an object or array block that follows a key in JSON-like text.
// Given src containing "\"KEY\": { ... }" or "\"KEY\": [ ... ]" and key "KEY",
// returns the substring between the matching braces/brackets (without the braces/brackets).
bool extract_object_block(const std::string& src, const std::string& key, std::string& out_block);

// Parse a number value for a given key inside src into out (double). Returns true on success.
bool parse_number(const std::string& src, const std::string& key, double& out);

// Parse a nested vec3 object under key into x, y, z doubles. Returns true on success.
bool parse_vec3(const std::string& src, const std::string& key, double& x, double& y, double& z);

// Parse a nested vec2i object under key with fields x and y as integers. Returns true on success.
bool parse_vec2i(const std::string& src, const std::string& key, int& outx, int& outy);

}

/** Lightweight 3D vector in world coordinates for rays. */
struct RayVec3 {
    double x;
    double y;
    double z;
};

/** Ray with position and direction in world coordinates. */
class Ray {
public:
	Ray();
	Ray(const RayVec3& position, const RayVec3& direction);

	/** Position/origin of the ray in world coordinates. */
	const RayVec3& getPosition() const;
	/** Direction of the ray in world coordinates (not normalized by default). */
	const RayVec3& getDirection() const;

	/** Set the position/origin of the ray in world coordinates. */
	void setPosition(const RayVec3& p);
	/** Set the direction of the ray in world coordinates. */
	void setDirection(const RayVec3& d);

private:
	RayVec3 m_position;
	RayVec3 m_direction;
};

// ---------------- Hit ----------------

/** Lightweight 3D vector in world coordinates (for hits). */
struct HitVec3 {
	double x;
	double y;
	double z;
};

/** Encapsulates information produced when a ray hits a surface. */
class Hit {
public:
	Hit(
		std::optional<HitVec3> intersectionPoint = std::nullopt,
		std::optional<double> distanceAlongRay = std::nullopt,
		std::optional<HitVec3> surfaceNormal = std::nullopt,
		std::optional<HitVec3> reflectedDirection = std::nullopt
	);

	// Getters (throw if unset)
	const HitVec3& getIntersectionPoint() const;
	double getDistanceAlongRay() const;
	const HitVec3& getSurfaceNormal() const;
	const HitVec3& getReflectedDirection() const;

	// Setters
	void setIntersectionPoint(const HitVec3& p);
	void setDistanceAlongRay(double t);
	void setSurfaceNormal(const HitVec3& n);
	void setReflectedDirection(const HitVec3& r);

	// Presence checks
	bool hasIntersectionPoint() const;
	bool hasDistanceAlongRay() const;
	bool hasSurfaceNormal() const;
	bool hasReflectedDirection() const;

private:
	std::optional<HitVec3> m_intersectionPoint;
	std::optional<double> m_distanceAlongRay;
	std::optional<HitVec3> m_surfaceNormal;
	std::optional<HitVec3> m_reflectedDirection;

	[[noreturn]] static void throwUnset(const char* name) {
		throw std::runtime_error(std::string("Hit: ") + name + " is unset");
	}
};

// ---------------- Image ----------------

/** Pixels are stored as 8-bit per channel; max value is preserved from input. */
struct Pixel {
	unsigned char r;
	unsigned char g;
	unsigned char b;
};

/** Image class for reading, manipulating, and writing PPM images. */
class Image {
public:
	explicit Image(const std::string& filename);
	Image(int width, int height);

	int width() const;
	int height() const;
	int maxValue() const;

	Pixel getPixel(int x, int y) const;
	void setPixel(int x, int y, int r, int g, int b);
	void addToPixel(int x, int y, int dr, int dg, int db);
	void setMaxValue(int maxVal);
	void write(const std::string& filename) const;

private:
	int m_width;
	int m_height;
	int m_maxValue;
	std::vector<Pixel> m_data;

	void loadFromFile(const std::string& filename);
};

#endif // UTILS_H


