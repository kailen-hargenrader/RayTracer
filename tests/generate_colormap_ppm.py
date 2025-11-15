import argparse
width, height = 480, 270
maxval = 255

def write_ppm(path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("P3\n")
        f.write(f"{width} {height}\n")
        f.write(f"{maxval}\n")
        for y in range(height):
            row = []
            for x in range(width):
                r = int(255 * x / (width - 1))
                g = int(255 * y / (height - 1))
                b = int(255 * ((x ^ y) % width) / (width - 1))
                row.append(f"{r} {g} {b}")
            f.write(" ".join(row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a 64x64 PPM (P3) test image.")
    parser.add_argument("--output", "-o", default="tests/raytracer/480x270_color_map.ppm", help="Output PPM file path")
    args = parser.parse_args()
    write_ppm(args.output)


if __name__ == "__main__":
    main()


