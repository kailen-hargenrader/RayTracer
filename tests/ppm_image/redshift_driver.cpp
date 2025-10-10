#include <iostream>
#include <string>
#include <exception>

#include "../../s2899221/Code/image.h"

// static inline unsigned char clamp_to_uchar(int v) {
//     if (v < 0) return 0;
//     if (v > 255) return 255;
//     return static_cast<unsigned char>(v);
// }

int main() {
    const std::string inputPath = "tests/ppm_image/reference_64x64.ppm";
    const std::string outputPath = "tests/ppm_image/redshift_reference_64x64.ppm";

    try {
        Image img(inputPath);

        // Redshift effect: boost red channel by a fixed delta, clamp at 255
        const int delta = 64;
        for (int y = 0; y < img.height(); ++y) {
            for (int x = 0; x < img.width(); ++x) {
                Pixel p = img.getPixel(x, y);
                const int r = static_cast<int>(p.r) + delta;
                img.setPixel(x, y, r, static_cast<int>(p.g), static_cast<int>(p.b));
            }
        }

        img.write(outputPath);
        std::cout << "Wrote redshifted image to: " << outputPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}


