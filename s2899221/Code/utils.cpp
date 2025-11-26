#include "utils.h"

#include <regex>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace util_json {

bool extract_object_block(const std::string& src, const std::string& key, std::string& out_block) {
	const std::string quoted_key = "\"" + key + "\"";
	const size_t key_pos = src.find(quoted_key);
	if (key_pos == std::string::npos) return false;

	size_t colon_pos = src.find(':', key_pos + quoted_key.size());
	if (colon_pos == std::string::npos) return false;

	// Find next non-space char after colon and determine opener
	size_t i = colon_pos + 1;
	while (i < src.size() && std::isspace(static_cast<unsigned char>(src[i]))) ++i;
	if (i >= src.size()) return false;
	char open_char = 0, close_char = 0;
	if (src[i] == '{' || src[i] == '[') {
		open_char = src[i];
		close_char = (open_char == '{') ? '}' : ']';
	} else {
		return false;
	}

	int depth = 0;
	bool in_string = false;
	bool is_escaped = false;
	size_t opener_pos = i;
	for (size_t j = opener_pos; j < src.size(); ++j) {
		const char c = src[j];
		if (in_string) {
			if (is_escaped) { is_escaped = false; }
			else if (c == '\\') { is_escaped = true; }
			else if (c == '"') { in_string = false; }
			continue;
		}
		if (c == '"') { in_string = true; continue; }
		if (c == open_char) { if (depth == 0) opener_pos = j; ++depth; }
		else if (c == close_char) {
			--depth;
			if (depth == 0) { out_block = src.substr(opener_pos + 1, j - opener_pos - 1); return true; }
		}
	}
	return false;
}

bool parse_number(const std::string& src, const std::string& key, double& out) {
	std::regex r("\\\"" + key + "\\\"\\s*:\\s*([-+eE0-9\\.]+)");
	std::smatch m; if (!std::regex_search(src, m, r)) return false; out = std::stod(m[1]); return true;
}

bool parse_vec3(const std::string& src, const std::string& key, double& x, double& y, double& z) {
	std::string block; if (!extract_object_block(src, key, block)) return false;
	double tx=0, ty=0, tz=0; if (!parse_number(block, "x", tx)) return false; if (!parse_number(block, "y", ty)) return false; if (!parse_number(block, "z", tz)) return false;
	x = tx; y = ty; z = tz; return true;
}

bool parse_vec2i(const std::string& src, const std::string& key, int& outx, int& outy) {
	std::string block; if (!extract_object_block(src, key, block)) return false;
	std::regex rx("\\\"x\\\"\\s*:\\s*([0-9]+)[^0-9]+\\\"y\\\"\\s*:\\s*([0-9]+)");
	std::smatch m; if (!std::regex_search(block, m, rx)) return false;
	outx = std::stoi(m[1]); outy = std::stoi(m[2]); return true;
}

}

// ---------------- Ray ----------------

Ray::Ray() : m_position{0,0,0}, m_direction{0,0,1} {}

Ray::Ray(const RayVec3& position, const RayVec3& direction)
	: m_position(position), m_direction(direction) {}

const RayVec3& Ray::getPosition() const { return m_position; }
const RayVec3& Ray::getDirection() const { return m_direction; }

void Ray::setPosition(const RayVec3& p) { m_position = p; }
void Ray::setDirection(const RayVec3& d) { m_direction = d; }

// ---------------- Hit ----------------

Hit::Hit(
	std::optional<HitVec3> intersectionPoint,
	std::optional<double> distanceAlongRay,
	std::optional<HitVec3> surfaceNormal,
	std::optional<HitVec3> reflectedDirection,
	std::optional<const Mesh*> mesh
)
	: m_intersectionPoint(intersectionPoint)
	, m_distanceAlongRay(distanceAlongRay)
	, m_surfaceNormal(surfaceNormal)
	, m_reflectedDirection(reflectedDirection)
	, m_mesh(mesh)
	, m_u(std::nullopt)
	, m_v(std::nullopt)
{}

const HitVec3& Hit::getIntersectionPoint() const {
	if (!m_intersectionPoint) throwUnset("intersectionPoint");
	return *m_intersectionPoint;
}

double Hit::getDistanceAlongRay() const {
	if (!m_distanceAlongRay) throwUnset("distanceAlongRay");
	return *m_distanceAlongRay;
}

const HitVec3& Hit::getSurfaceNormal() const {
	if (!m_surfaceNormal) throwUnset("surfaceNormal");
	return *m_surfaceNormal;
}

const HitVec3& Hit::getReflectedDirection() const {
	if (!m_reflectedDirection) throwUnset("reflectedDirection");
	return *m_reflectedDirection;
}

const Mesh* Hit::getMesh() const {
	if (!m_mesh) throwUnset("mesh");
	return *m_mesh;
}

void Hit::setIntersectionPoint(const HitVec3& p) { m_intersectionPoint = p; }
void Hit::setDistanceAlongRay(double t) { m_distanceAlongRay = t; }
void Hit::setSurfaceNormal(const HitVec3& n) { m_surfaceNormal = n; }
void Hit::setReflectedDirection(const HitVec3& r) { m_reflectedDirection = r; }
void Hit::setMesh(const Mesh* m) { m_mesh = m; }

bool Hit::hasIntersectionPoint() const { return m_intersectionPoint.has_value(); }
bool Hit::hasDistanceAlongRay() const { return m_distanceAlongRay.has_value(); }
bool Hit::hasSurfaceNormal() const { return m_surfaceNormal.has_value(); }
bool Hit::hasReflectedDirection() const { return m_reflectedDirection.has_value(); }
bool Hit::hasMesh() const { return m_mesh.has_value(); }

void Hit::setUV(double u, double v) { m_u = u; m_v = v; }
bool Hit::hasUV() const { return m_u.has_value() && m_v.has_value(); }
double Hit::getU() const { if (!m_u) throwUnset("u"); return *m_u; }
double Hit::getV() const { if (!m_v) throwUnset("v"); return *m_v; }
// ---------------- Image ----------------

static void skipComments(std::istream& in) {
	int c = in.peek();
	while (std::isspace(c) || c == '#') {
		if (c == '#') { std::string dummy; std::getline(in, dummy); }
		else { in.get(); }
		c = in.peek();
	}
}

Image::Image(const std::string& filename)
	: m_width(0), m_height(0), m_maxValue(255) { loadFromFile(filename); }

Image::Image(int width, int height)
	: m_width(width), m_height(height), m_maxValue(255), m_data(static_cast<size_t>(width) * static_cast<size_t>(height)) {}

int Image::width() const { return m_width; }
int Image::height() const { return m_height; }
int Image::maxValue() const { return m_maxValue; }

Pixel Image::getPixel(int x, int y) const {
	return m_data[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)];
}

static inline unsigned char clamp_to_max(int v, int maxv) {
	if (v < 0) return 0;
	if (v > maxv) return static_cast<unsigned char>(maxv);
	return static_cast<unsigned char>(v);
}

void Image::setPixel(int x, int y, int r, int g, int b) {
	Pixel p; const int maxv = m_maxValue;
	p.r = clamp_to_max(r, maxv);
	p.g = clamp_to_max(g, maxv);
	p.b = clamp_to_max(b, maxv);
	m_data[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] = p;
}

void Image::addToPixel(int x, int y, int dr, int dg, int db) {
	Pixel p = getPixel(x, y); const int maxv = m_maxValue;
	p.r = clamp_to_max(static_cast<int>(p.r) + dr, maxv);
	p.g = clamp_to_max(static_cast<int>(p.g) + dg, maxv);
	p.b = clamp_to_max(static_cast<int>(p.b) + db, maxv);
	m_data[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] = p;
}

void Image::setMaxValue(int maxVal) {
	if (maxVal < 1 || maxVal > 255) throw std::runtime_error("Invalid PPM max value, must be [1,255].");
	m_maxValue = maxVal;
}

void Image::write(const std::string& filename) const {
	std::ofstream out(filename, std::ios::binary);
	if (!out) { throw std::runtime_error("Failed to open file for writing: " + filename); }
	out << "P6\n" << m_width << " " << m_height << "\n" << m_maxValue << "\n";
	std::vector<unsigned char> buf;
	buf.reserve(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3u);
	for (const Pixel& p : m_data) {
		buf.push_back(clamp_to_max(p.r, m_maxValue));
		buf.push_back(clamp_to_max(p.g, m_maxValue));
		buf.push_back(clamp_to_max(p.b, m_maxValue));
	}
	out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
	if (!out) { throw std::runtime_error("Failed while writing image data to: " + filename); }
}

void Image::loadFromFile(const std::string& filename) {
	std::ifstream in(filename, std::ios::binary);
	if (!in) { throw std::runtime_error("Failed to open image file: " + filename); }

	std::string magic; in >> magic;
	if (magic != "P6" && magic != "P3") { throw std::runtime_error("Unsupported PPM format (expected P3 or P6): " + magic); }

	skipComments(in);
	in >> m_width; skipComments(in);
	in >> m_height; skipComments(in);
	int maxval = 0; in >> maxval; in.get();
	if (m_width <= 0 || m_height <= 0 || maxval <= 0) { throw std::runtime_error("Invalid PPM header values"); }
	if (maxval > 255) { throw std::runtime_error("Unsupported PPM max value (>255 not supported)"); }

	m_maxValue = maxval;
	m_data.resize(static_cast<size_t>(m_width) * static_cast<size_t>(m_height));

	if (magic == "P6") {
		const size_t numBytes = static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3u;
		std::vector<unsigned char> buf(numBytes);
		in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
		if (static_cast<size_t>(in.gcount()) != numBytes) { throw std::runtime_error("Unexpected EOF while reading binary PPM data"); }
		for (size_t i = 0, p = 0; i < m_data.size(); ++i) { m_data[i].r = buf[p++]; m_data[i].g = buf[p++]; m_data[i].b = buf[p++]; }
	} else {
		for (size_t i = 0; i < m_data.size(); ++i) {
			int r, g, b;
			skipComments(in);
			if (!(in >> r)) throw std::runtime_error("Unexpected EOF reading P3 data (r)");
			skipComments(in);
			if (!(in >> g)) throw std::runtime_error("Unexpected EOF reading P3 data (g)");
			skipComments(in);
			if (!(in >> b)) throw std::runtime_error("Unexpected EOF reading P3 data (b)");
			if (r < 0 || r > maxval || g < 0 || g > maxval || b < 0 || b > maxval) { throw std::runtime_error("P3 pixel out of range"); }
			m_data[i].r = static_cast<unsigned char>(r);
			m_data[i].g = static_cast<unsigned char>(g);
			m_data[i].b = static_cast<unsigned char>(b);
		}
	}
}

