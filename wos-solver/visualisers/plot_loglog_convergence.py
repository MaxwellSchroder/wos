import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Load CSV
data = np.loadtxt("../error_plot_x0.csv", delimiter=',')
n_walks = data[:, 0]
l1_error = data[:, 1]

# Compute logs
log_n = np.log10(n_walks)
log_l1 = np.log10(l1_error)

# Perform linear regression: log(L1 error) = slope * log(nWalks) + intercept
slope, intercept, r_value, p_value, std_err = linregress(log_n, log_l1)
fitted_line = slope * log_n + intercept

# Calculate reference line with slope -0.5 passing through first data point
ref_slope = -0.5
ref_intercept = log_l1[0] - ref_slope * log_n[0]
ref_line = ref_slope * log_n + ref_intercept

plt.figure(figsize=(8, 5))

# Plot actual data and fitted line
plt.plot(log_n, log_l1, 'o-', label='Data (log-log)', markersize=4)
plt.plot(log_n, fitted_line, 'r--', label=f'Fit: slope = {slope:.3f}')

# Plot reference convergence line
plt.plot(log_n, ref_line, 'k--', label='Reference: slope = -0.5', alpha=0.6)

plt.xlabel("log₁₀(Number of Walks)", fontsize=12)
plt.ylabel("log₁₀(L1 Error)", fontsize=12)
plt.title("Log-Log Convergence of WoS at Single Point", fontsize=14)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("loglog_convergence_plot.png", dpi=300)
plt.show()

# Print convergence rate
print(f"Estimated convergence rate (slope): {slope:.4f}")
print(f"R² = {r_value**2:.4f}")
