import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy.stats import linregress

# Load CSV, skip header row
data = np.loadtxt("../all_epsilon_convergence.csv", delimiter=',', skiprows=1)

epsilons = data[:, 0]
n_walks = data[:, 1]
l1_errors = data[:, 2]
cumulative_time = data[:,3]

# Find unique epsilon values
unique_epsilons = np.unique(epsilons)[::-1]  # Reverse the order
print(str(unique_epsilons))

# Create the plot
# plt.figure(figsize=(8, 5))

def show_time_vs_epsilon():
    # Create a dictionary to track max n_walk per epsilon
    epsilon_to_stats = {}

    # Track max n_walk per epsilon, and save time and error
    epsilon_to_stats = {}

    for i in range(len(data)):
        eps = epsilons[i]
        n = n_walks[i]
        if eps not in epsilon_to_stats or n > epsilon_to_stats[eps][1]:
            epsilon_to_stats[eps] = (cumulative_time[i], n, l1_errors[i])

    # Extract sorted values
    sorted_epsilons = sorted(epsilon_to_stats.keys(), reverse=True)
    execution_times = [epsilon_to_stats[eps][0] for eps in sorted_epsilons]
    l1_errors_sorted = [epsilon_to_stats[eps][2] for eps in sorted_epsilons]
    
    # Linear Regression on values
    # Use log-log scale
    log_eps = np.log(sorted_epsilons)
    log_times = np.log(execution_times)
    
    # Perform linear regression
    slope, intercept, r_value, _, _ = linregress(log_eps, log_times)

    # Compute predicted line
    log_pred = slope * log_eps + intercept
    predicted_times = np.exp(log_pred)
    
    fig, ax1 = plt.subplots()
    
    # Execution time plot
    ax1.plot(sorted_epsilons, execution_times, 'bo-', label='Execution Time')
    ax1.plot(sorted_epsilons, predicted_times, 'b--', label=f"Fit: time ≈ {np.exp(intercept):.2f} * ε^{slope:.2f}")
    ax1.set_xscale('log')
    ax1.invert_xaxis()
    ax1.set_xlabel("Epsilon shell (log scale, decreasing)")
    ax1.set_ylabel("Execution Time (s)", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.set_xticks(sorted_epsilons)  # Explicitly set ticks
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
    
    # L1 error on second axis
    # Convert l1_errors to log scale for ticks
    ax2 = ax1.twinx()
    ax2.plot(sorted_epsilons, l1_errors_sorted, 'ro-', label='L1 Error')
    ax2.set_ylabel("Log L1 Error", color='r')
    ax2.set_yscale('log')
    min_val = min(l1_errors_sorted)
    max_val = max(l1_errors_sorted)
    
    # Create 6 log-spaced ticks between min and max 
    log_ticks = np.logspace(np.log10(min_val), np.log10(max_val), num=6)
    ax2.set_yticks(log_ticks)
    ax2.set_yticklabels([f"$10^{{{round(val,2)}}}$" for val in log_ticks])
    ax2.tick_params(axis='y', labelcolor='r')
    
    # Combine legends
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc='upper left')
    
    ## Display
    plt.title("For nWalks = 10^2: Total Execution Time and L1 Error vs Epsilon Shell Size")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# show_time_vs_epsilon()

def show_L1_error_vs_time():
    fig, ax1 = plt.subplots()
    
    for eps in unique_epsilons:
        mask = epsilons == eps
        ax1.plot(
            cumulative_time[mask],
            l1_errors[mask],
            marker='o',
            linestyle='-',
            linewidth=1,
            markersize=2,
            label=f"ε = {eps:.3g}"
        )
        
    # X-Axis
    plt.xlabel("Time (seconds)", fontsize=12)
    
    # Y-Axis
    plt.yscale('log')
    plt.ylabel("Log L1 Error (N)", fontsize=12)
    yUpperLim= 10**2
    yLowerlim = 10**-2
    plt.ylim([yLowerlim,yUpperLim])
    
    ## Display
    plt.legend(title="Epsilon Values", fontsize=10, title_fontsize=11)
    plt.title("L1 Error vs Execution Time")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()
    
show_L1_error_vs_time()


# for eps in unique_epsilons:
#     mask = epsilons == eps
#     plt.plot(
#         n_walks[mask],
#         l1_errors[mask],
#         marker='o',
#         linestyle='-',
#         linewidth=1,
#         markersize=2,
#         label=f"ε = {eps:.3g}"
#     )

# for eps in unique_epsilons:
#     mask = epsilons == eps
#     plt.plot(
#         n_walks[mask],
#         cumulative_time[mask],
#         marker='o',
#         linestyle='-',
#         linewidth=1,
#         markersize=2,
#         label=f"ε = {eps:.3g}"
#     )

# Label axes
# plt.xlabel("Number of Walks (N)", fontsize=12)
# plt.ylabel("Time (seconds)", fontsize=12)
# plt.title("Plot of Execution Time vs Number of Walks at a Single Point\nfor Different ε Values", fontsize=14)


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

# plt.ylim([0, 20])  # adjust if needed
# plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Add legend
# plt.legend(title="Epsilon Values", fontsize=10, title_fontsize=11)

# plt.tight_layout()
# plt.savefig("l1_convergence_multiple_epsilons.png", dpi=300)
# plt.savefig("Log_Log_l1_convergence_multiple_epsilons.png", dpi=300)
# plt.show()
