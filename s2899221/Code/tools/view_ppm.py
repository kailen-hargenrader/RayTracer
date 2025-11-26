import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python view_ppm.py <image.ppm>")
        sys.exit(1)
    path = sys.argv[1]
    try:
        from PIL import Image  # pip install pillow
    except ImportError:
        print("Pillow not installed. Install with: pip install pillow")
        sys.exit(1)
    img = Image.open(path)
    img.show()

if __name__ == "__main__":
    main()


