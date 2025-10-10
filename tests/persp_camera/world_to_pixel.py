import numpy as np

# Camera parameters
camera_location = np.array([0.0, -5.0, 2.0])
camera_direction = np.array([0.0, 0.9578263213120298, -0.2873477652633898])
focal_length_mm = 50.0
sensor_width_mm = 36.0
sensor_height_mm = 24.0
image_width_px = 1920
image_height_px = 1080

# World point to project
point_world = np.array([1.0, 2.0, 0.5])

# Normalize the forward (camera direction) vector
forward = camera_direction / np.linalg.norm(camera_direction)

# Assume world up is Z-up
world_up = np.array([0.0, 0.0, 1.0])

# Compute right and up vectors
right = np.cross(world_up, forward)
right = right / np.linalg.norm(right)

# Handle edge case: camera looking directly up/down
if np.linalg.norm(right) < 1e-6:
    world_up = np.array([0.0, 1.0, 0.0])
    right = np.cross(world_up, forward)
    right = right / np.linalg.norm(right)

up = np.cross(forward, right)
up = up / np.linalg.norm(up)

# Build rotation matrix (camera-to-world basis)
R = np.column_stack((right, up, forward))

# Transform world point into camera space
vec = point_world - camera_location
camera_coords = R.T @ vec
Xc, Yc, Zc = camera_coords

# Project onto sensor plane
x_sensor = focal_length_mm * (Xc / Zc)
y_sensor = focal_length_mm * (Yc / Zc)

# Convert to pixel coordinates
x_pixel = ((x_sensor / sensor_width_mm) + 0.5) * image_width_px
y_pixel = (1.0 - ((y_sensor / sensor_height_mm) + 0.5)) * image_height_px

print(f"Pixel coordinates: ({int(round(x_pixel))}, {int(round(y_pixel))})")
