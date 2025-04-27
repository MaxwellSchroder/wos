import matplotlib.pyplot as plt
import numpy as np

# Load CSV, skip header row
data = np.loadtxt("../all_epsilon_convergence.csv", delimiter=',', skiprows=1)

epsilons = data[:, 0]
n_walks = data[:, 1]
l1_errors = data[:, 2]

# Find unique epsilon values
unique_epsilons = np.unique(epsilons)[::-1]  # Reverse the order
print(str(unique_epsilons))

# Create the plot
plt.figure(figsize=(8, 5))

for eps in unique_epsilons:
    mask = epsilons == eps
    plt.plot(
        n_walks[mask],
        l1_errors[mask],
        marker='o',
        linestyle='-',
        linewidth=1,
        markersize=2,
        label=f"ε = {eps:.3g}"
    )

# Label axes
plt.xlabel("Number of Walks (N)", fontsize=12)
plt.ylabel("L1 Error", fontsize=12)
plt.title("Plot of Convergence of Walk on Spheres at a Single Point\nfor Different ε Values", fontsize=14)


# # Log Label Axes
# plt.xlabel("Log Number of Walks (N)", fontsize=12)
# plt.ylabel("Log L1 Error", fontsize=12)
# plt.title("Log-Log Plot of Convergence of Walk on Spheres at a Single Point\nfor Different ε Values", fontsize=14)
# # Optional: log-log plot
# plt.xscale('log')
# plt.yscale('log')
# # Reference convergence line (O(1/√N))
# ref_x = np.array([min(n_walks), max(n_walks)])
# ref_y = l1_errors[0] * (ref_x / n_walks[0])**(-0.5)
# plt.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

plt.ylim([0, 20])  # adjust if needed
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Add legend
plt.legend(title="Epsilon Values", fontsize=10, title_fontsize=11)

plt.tight_layout()
plt.savefig("l1_convergence_multiple_epsilons.png", dpi=300)
plt.show()
