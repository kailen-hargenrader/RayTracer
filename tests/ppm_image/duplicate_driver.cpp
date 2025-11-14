#include <iostream>
#include <string>
#include <exception>

#include "../../s2899221/Code/utils.h"

int main() {
    const std::string inputPath = "tests/ppm_image/reference_64x64.ppm";
    const std::string outputPath = "tests/ppm_image/duplicate_reference_64x64.ppm";

    try {
        Image img(inputPath);
        img.write(outputPath);
        std::cout << "Wrote duplicate image to: " << outputPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}


