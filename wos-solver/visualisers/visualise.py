import numpy as np
import matplotlib.pyplot as plt
import sys


data = np.loadtxt('out.csv', delimiter=',')
imageSize = 256

plt.figure(figsize=(8, 6))
plt.imshow(data, cmap='inferno', interpolation='nearest', vmin=0)
plt.colorbar(label='Temperature (Kelvin)')
plt.title("\n Heat Distribution on 2D Shape. \n Visualisation of WoS Algorithm \n N walks = 128 ")
plt.axis([0,imageSize,0,imageSize])
plt.xlabel("X Position")

# Set custom ticks and labels for the x-axis and y-axis
# The ticks are scaled to match the imageSize-pixel range
pixel_positions = [0.2 * imageSize, 0.4 * imageSize, 0.6 * imageSize, 0.8 * imageSize, 1.0 * imageSize]
scale_labels = [0.2, 0.4, 0.6, 0.8, 1.0]  # The labels for the ticks

# Set the ticks and labels for both axes
plt.xticks(pixel_positions, scale_labels)
plt.yticks(pixel_positions, scale_labels)

plt.ylabel("Y Position")
plt.show()

