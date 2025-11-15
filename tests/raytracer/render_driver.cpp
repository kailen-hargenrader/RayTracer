// render_driver.cpp
#include <iostream>
#include <string>

#include "../../s2899221/Code/raytracer.h"
// Pull implementations directly so this driver links as a single TU
#include "../../s2899221/Code/utils.cpp"
#include "../../s2899221/Code/mesh.cpp"
#include "../../s2899221/Code/persp_camera.cpp"
#include "../../s2899221/Code/raytracer.cpp"

static void print_usage() {
	std::cerr << "Usage: render_driver <scene.json> <camera_id> <output.ppm>\n";
	std::cerr << "- scene.json: path to JSON scene file\n";
	std::cerr << "- camera_id: id string of camera to render\n";
	std::cerr << "- output.ppm: destination PPM file path\n";
}

int main(int argc, char** argv) {
	if (argc < 4 || argc > 5) {
		print_usage();
		return 1;
	}

	const std::string json_path = argv[1];
	const std::string camera_id = argv[2];
	const std::string output_path = argv[3];
	const std::string color_ppm_path = argv[4];
	

	RayTracer rt;
	if (!rt.load_from_json(json_path)) {
		std::cerr << "Failed to load scene: " << json_path << "\n";
		return 1;
	}

	if (!rt.render_unaccelerated_ppm(camera_id, output_path, color_ppm_path)) {
		std::cerr << "Render failed (bad camera id or write error)\n";
		return 1;
	}

	std::cout << "Wrote image to: " << output_path << "\n";
	return 0;
}

