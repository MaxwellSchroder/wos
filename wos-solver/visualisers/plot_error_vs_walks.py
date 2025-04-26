import matplotlib.pyplot as plt
import numpy as np

# Load the CSV data
data = np.loadtxt("../error_plot_x0.csv", delimiter=',')
n_walks = data[:, 0]
l1_error = data[:, 1]

# Plot the data
plt.figure(figsize=(8, 5))
plt.plot(n_walks, l1_error, marker='o', linestyle='none', linewidth=1.5)
plt.ylim(0, max(l1_error) * 1.1)
# plt.plot(n_walks, l1_error, marker='o', linestyle='-', linewidth=1.5)

# Label axes
plt.xlabel("Number of Walks (N)", fontsize=12)
plt.ylabel("L1 Error", fontsize=12)
plt.title("Convergence of Walk on Spheres at a Single Point", fontsize=14)

# Optional: Uncomment to use log-log scale (common in Monte Carlo convergence plots)
# plt.xscale('log')
# plt.yscale('log')

plt.grid(True, which='both', linestyle='--', linewidth=0.2)
plt.tight_layout()
plt.savefig("l1_convergence_plot.png", dpi=300)
plt.show()
