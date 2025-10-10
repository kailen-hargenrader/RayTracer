#ifndef IMAGE_H
#define IMAGE_H

#include <string>
#include <vector>

/**
 * Pixels are stored as 8-bit per channel and use the same scale as the
 * image's max value. This implementation supports PPM max color values in
 * the range [1, 255].
 */
struct Pixel {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

/**
 * Image class for reading, manipulating, and writing PPM images.
 *
 * Supports P3 (ASCII) and P6 (binary) input. Output is P6. The stored
 * max color value is preserved from input (or defaults to 255) and is
 * used for clamping and when writing.
 */
class Image {
public:
    /** Construct by loading a PPM file from disk (P3 or P6). Throws on failure. */
    explicit Image(const std::string& filename);

    /** Construct an empty image with given dimensions, filled with black. */
    Image(int width, int height);

    int width() const;
    int height() const;
    /** Returns the PPM max color value associated with this image (1..255). */
    int maxValue() const;

    /** Get a pixel at (x, y). No bounds checking. */
    Pixel getPixel(int x, int y) const;

    /** Set a pixel using integer RGB with clamping to [0,maxValue()]. */
    void setPixel(int x, int y, int r, int g, int b);

    /** Add deltas to a pixel with clamping to [0,maxValue()]. Useful for effects. */
    void addToPixel(int x, int y, int dr, int dg, int db);

    /**
     * Set the max color value used for clamping and when writing. Must be
     * in [1, 255]. Values outside this range are invalid.
     */
    void setMaxValue(int maxVal);

    /** Write to a PPM file (P6 binary). Throws on failure. Uses current max value. */
    void write(const std::string& filename) const;

private:
    int m_width;
    int m_height;
    int m_maxValue; // PPM max color value (1..255). 16-bit PPM is not supported.
    std::vector<Pixel> m_data;

    // Helpers
    void loadFromFile(const std::string& filename);
};

#endif // IMAGE_H


