from PIL import Image
from pathlib import Path
import sys


def main() -> None:
    # Allow optional CLI arg: path to PPM; default to file next to this script
    if len(sys.argv) > 1:
        img_path = Path(sys.argv[1])
    else:
        img_path = Path(__file__).parent / "rendered_cube_cylinder_plane_sphere.ppm"

    if not img_path.exists():
        print(f"PPM not found: {img_path}")
        sys.exit(1)

    img = Image.open(img_path)
    img.show()


if __name__ == "__main__":
    main()