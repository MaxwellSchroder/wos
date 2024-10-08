import numpy as np
import matplotlib.pyplot as plt

# Number of points on the lines and curves
n_points_flat = 50  # Number of points on each flat line
n_points_curve = 100  # Number of points on each curved line

# 1. Generate points for the bottom flat line (from (0.2,0.2) to (0.4,0.2))
x_bottom = np.linspace(0.2, 0.4, n_points_flat)
y_bottom = np.full_like(x_bottom, 0.2)  # y = 0.2 for the entire flat bottom line
temperature_bottom = np.full_like(x_bottom, fill_value=150) # Create a temperature field, set to 150K

# 2. Generate points for the right flat line (from (0.8,0.6) to (0.8,0.8))
x_right = np.full(n_points_flat, 0.8)  # x = 0.8 for the entire right line
y_right = np.linspace(0.6, 0.8, n_points_flat)
temperature_right = np.full_like(x_right, fill_value=150) # Create a temperature field, set to 150K


# 3. Generate points for the curved lines (circular arc or parametric curve)
# First curved line connecting the bottom flat line to the start of the right flat line
theta1 = np.linspace(0, np.pi / 2, n_points_curve)  # Angle from 0 to 90 degrees

# SMALLER curve parameters (using an arbitrary circular arc for the bend)
radius = 0.4  # You can adjust this radius for the curve
x_curve1 = 0.8 - radius * np.cos(theta1)
y_curve1 = 0.2 + radius * np.sin(theta1)
temperature_small_curve = np.full_like(theta1, fill_value=100) # Create a temperature field, set to 100K



# Larger curve - has to be reversed as it needs to come down for the points to continue to be anti-clockwise winding order
large_radius = 0.8 - 0.2
x_curve2 = 0.8 - large_radius * np.cos(theta1)
y_curve2 = 0.2 + large_radius * np.sin(theta1)
x_curve2 = x_curve2[::-1]
y_curve2 = y_curve2[::-1]
temperature_larger_curve = np.full_like(theta1, fill_value=300) # Create a temperature field, set to 300K
# reverse temperature for later
temperature_larger_curve = temperature_larger_curve[::-1]


# combine all arrays
x_combined = np.concatenate([x_bottom, x_curve1, x_right, x_curve2])
y_combined = np.concatenate([y_bottom, y_curve1, y_right, y_curve2])
temperature_combined = np.concatenate([temperature_bottom, temperature_small_curve, temperature_right, temperature_larger_curve])


# Plot the geometry
# plt.plot(x_bottom, y_bottom, 'b-', label='Bottom Line')
# plt.plot(x_right, y_right, 'b-', label='Right Line')
# plt.plot(x_curve1, y_curve1, 'r-', label="Bottom curve")
# plt.plot(x_curve2, y_curve2, 'r-', label="Larger curve")
plt.plot(x_combined, y_combined, 'r-')


plt.axis([0,1,0,1])
plt.gca().set_aspect('equal', adjustable='box')  # This ensures equal scaling
# plt.scatter(x_points, y_points, color='red', s=5)  # Show individual points
plt.title('2D Geometry with Flat and Curved Lines')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.grid(True)
plt.legend()
# plt.show()


# Write to a file where each line contains x and y comma-separated
with open('combined_coordinates.csv', 'w') as f:
    for x, y, t in zip(x_combined, y_combined, temperature_combined):
        f.write(f"{x},{y},{t}\n")

print("Coordinates have been written to combined_coordinates.csv.")

# print(str(temperature_combined))

