import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy.stats import linregress

# Load CSV, skip header row
data = np.loadtxt("../full_solves/all_epsilon_convergence_7.csv", delimiter=',', skiprows=1)
# data = np.loadtxt("../all_epsilon_convergence.csv", delimiter=',', skiprows=1)

# data = np.loadtxt("../interior_point_experiment/dirichlet_solve2.csv", delimiter=',', skiprows=1)
# data = np.loadtxt("../interior_point_experiment/centroid_solve.csv", delimiter=',', skiprows=1)
# data = np.loadtxt("../interior_point_experiment/neumann_solve2.csv", delimiter=',', skiprows=1)

epsilons = data[:, 0]
n_walks = data[:, 1]
l1_errors = data[:, 2]
cumulative_time = data[:,3]
steps = data[:,4]

hardcoded_true_T_solution = 366.6345156816051

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

def show_percent_error_vs_time():
    hardcoded_true_T_solution = 366.6345156816051

    fig, ax1 = plt.subplots()

    for eps in unique_epsilons:
        mask = epsilons == eps
        percent_error = (l1_errors[mask] / hardcoded_true_T_solution) * 100

        ax1.plot(
            cumulative_time[mask],
            percent_error,
            marker='o',
            linestyle='-',
            linewidth=1,
            markersize=2,
            label=f"ε = {eps:.3g}"
        )

    # X-Axis
    ax1.set_xlabel("Time (seconds)", fontsize=12)

    # Y-Axis: linear % error
    ax1.set_ylabel("Relative Error (%)", fontsize=12)
    ax1.set_ylim(0, min(1, np.max(percent_error) * 1.1))  # auto-scale but clamp below 10% range

    # Title, legend, grid
    ax1.set_title("Percent Error vs Execution Time")
    ax1.legend(title="Epsilon Values", fontsize=10, title_fontsize=11)
    ax1.grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.show()

def show_l1_and_percent_error_vs_time():
    hardcoded_true_T_solution = 366.6345156816051

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1]})

    # --- Upper plot: Log L1 Error ---
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

    ax1.set_yscale('log')
    ax1.set_ylabel("Log L1 Error", fontsize=12)
    ax1.set_title("Log L1 Error vs Execution Time")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax1.legend(title="Epsilon", fontsize=10, title_fontsize=11)
    ax1.tick_params(labelbottom=True)

    # --- Lower plot: Percent Relative Error ---
    for eps in unique_epsilons:
        mask = epsilons == eps
        percent_error = (l1_errors[mask] / hardcoded_true_T_solution) * 100
        ax2.plot(
            cumulative_time[mask],
            percent_error,
            marker='o',
            linestyle='-',
            linewidth=1,
            markersize=2,
            label=f"ε = {eps:.3g}"
        )

    ax2.set_ylabel("Relative Error (%)", fontsize=12)
    ax2.set_xlabel("Time (seconds)", fontsize=12)
    ax2.set_ylim(0, min(1, np.max(percent_error) * 1.1))
    ax2.set_title("Percent Relative Error vs Execution Time")
    ax2.grid(True, linestyle='--', linewidth=0.5)

    # --- Layout and caption ---
    fig.suptitle("Comparison of Absolute and Relative Error Over Time", fontsize=14)
    fig.text(
        0.5, 0.01,
        "Top panel shows absolute L1 error on a log scale.\nBottom panel shows relative error as a percentage of the true value (T = 366.63 K).",
        ha='center', fontsize=10
    )

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, bottom=0.18)
    plt.show()

def show_L1_error_vs_N_Walks():
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 10), sharex=True)
    
    for eps in unique_epsilons:
        mask = epsilons == eps
        
        x = n_walks[mask]
        y = l1_errors[mask]
        
        # Plot on both subplots
        ax1.plot(x, y, marker='o', linestyle='-', linewidth=1, markersize=2, label=f"ε = {eps:.3g}")
        ax2.plot(x, y, marker='o', linestyle='-', linewidth=1, markersize=2, label=f"ε = {eps:.3g}")
    
    ## Log-Log Plot (Top)
    ax1.set_ylabel("Log L1 Error")
    ax1.set_xlabel("Number of Walks (N)")
    ax1.tick_params(labelbottom=True)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_title("Log-Log Convergence of WoSt at a Single Point")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    # Reference convergence line (O(1/√N))
    ref_x = np.array([min(n_walks), max(n_walks)])
    ref_y = l1_errors[0] * (ref_x / n_walks[0])**(-0.5)
    ax1.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

    ## Linear Y-Axis Plot (Bottom)
    ax2.set_xlabel("Number of Walks (N)")
    ax2.set_ylabel("L1 Error")
    ax2.set_xscale('log')
    ax2.set_yscale('linear')
    ax2.set_ylim([-1,2])
    ax2.set_title("Log-X Linear-Y Plot of Convergence (Same Data)")
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

    ## Display
    ax1.legend(title="Epsilon Values", fontsize=9, title_fontsize=10)
    # plt.title("Log-Log Plot of Convergence of Walk on Spheres\n at a Single Point for Different ε Values", fontsize=14)
    # plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    ## SAVE
    # plt.savefig("l1_convergence_multiple_epsilons.png", dpi=300)
    # plt.savefig("Log_Log_l1_convergence_multiple_epsilons.png", dpi=300)

def show_L1_error_vs_N_Walks_single():
    fig, ax1 = plt.subplots(figsize=(7, 4))

    for eps in unique_epsilons:
        mask = epsilons == eps

        x = n_walks[mask]
        y = l1_errors[mask]

        ax1.plot(x, y, marker='o', linestyle='-', linewidth=1, markersize=2, label=f"ε = {eps:.3g}")

    # Log-Log Plot
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Log Number of Walks (N)")
    ax1.set_ylabel("Log L1 Error")
    ax1.set_title("Log-Log Convergence of WoS at a Single Point")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Reference convergence line (O(1/√N))
    ref_x = np.array([min(n_walks), max(n_walks)])
    ref_y = l1_errors[0] * (ref_x / n_walks[0])**(-0.5)
    ax1.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

    ax1.legend(title="Epsilon Values", fontsize=9, title_fontsize=10)

    plt.tight_layout()
    plt.show()

def show_L1_error_vs_N_Walks_multiple_points():
    # Load datasets for 3 interior points
    centroid_data = np.loadtxt("../interior_point_experiment/centroid_solve.csv", delimiter=',', skiprows=1)
    dirichlet_data = np.loadtxt("../interior_point_experiment/dirichlet_solve.csv", delimiter=',', skiprows=1)
    neumann_data = np.loadtxt("../interior_point_experiment/neumann_solve.csv", delimiter=',', skiprows=1)

    # Unpack columns: assuming consistent structure [epsilon, n_walks, l1_error, time, steps, ...]
    centroid_n = centroid_data[:, 1]
    centroid_l1 = centroid_data[:, 2]

    dirichlet_n = dirichlet_data[:, 1]
    dirichlet_l1 = dirichlet_data[:, 2]

    neumann_n = neumann_data[:, 1]
    neumann_l1 = neumann_data[:, 2]

    # Plot
    fig, ax1 = plt.subplots(figsize=(7, 4))

    ax1.plot(centroid_n, centroid_l1, 'o-', linewidth=1, markersize=2, label="Centroid Point")
    ax1.plot(dirichlet_n, dirichlet_l1, 's-', linewidth=1, markersize=2, label="Near Dirichlet Boundary")
    ax1.plot(neumann_n, neumann_l1, '^-', linewidth=1, markersize=2, label="Near Neumann Boundary")

    # Log-Log Axes
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Number of Walks (N)")
    ax1.set_ylabel("Log L1 Error")
    ax1.set_title("Log-Log Convergence at Different Interior Points. Eps = 0.00125")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Reference convergence line
    ref_x = np.array([min(centroid_n), max(centroid_n)])
    ref_y = centroid_l1[0] * (ref_x / centroid_n[0])**(-0.5)
    ax1.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

    ax1.legend(title="Interior Evaluation Points", fontsize=9, title_fontsize=10)
    plt.tight_layout()
    plt.show()
    
def show_L1_error_vs_N_Walks_multiple_seeds():
    # Load datasets for 3 interior points
    seed1234_data = np.loadtxt("../seed/seed1234.csv", delimiter=',', skiprows=1)
    seed5678_data = np.loadtxt("../seed/seed5678.csv", delimiter=',', skiprows=1)
    seed9123_data = np.loadtxt("../seed/seed9123.csv", delimiter=',', skiprows=1)
    seed4567_data = np.loadtxt("../seed/seed4567.csv", delimiter=',', skiprows=1)
    seed8912_data = np.loadtxt("../seed/seed8912.csv", delimiter=',', skiprows=1)

    # Unpack columns: assuming consistent structure [epsilon, n_walks, l1_error, time, steps, ...]
    seed1234_n = seed1234_data[:, 1]
    seed1234_l1 = seed1234_data[:, 2]

    seed5678_n = seed5678_data[:, 1]
    seed5678_l1 = seed5678_data[:, 2]
    
    seed9123_n = seed9123_data[:, 1]
    seed9123_l1 = seed9123_data[:, 2]
    
    seed4567_n = seed4567_data[:, 1]
    seed4567_l1 = seed4567_data[:, 2]
    
    seed8912_n = seed8912_data[:, 1]
    seed8912_l1 = seed8912_data[:, 2]

    # Plot
    fig, ax1 = plt.subplots(figsize=(7, 4))

    ax1.plot(seed1234_n, seed1234_l1, 'o-', linewidth=1, markersize=2, label="Seed 1234")
    ax1.plot(seed5678_n, seed5678_l1, 's-', linewidth=1, markersize=2, label="Seed 5678")
    ax1.plot(seed9123_n, seed9123_l1, 's-', linewidth=1, markersize=2, label="Seed 9123")
    ax1.plot(seed4567_n, seed4567_l1, 's-', linewidth=1, markersize=2, label="Seed 4567")
    ax1.plot(seed8912_n, seed8912_l1, 's-', linewidth=1, markersize=2, label="Seed 8912")

    # Log-Log Axes
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Number of Walks (N)")
    ax1.set_ylabel("Log L1 Error")
    ax1.set_title("Log-Log Convergence at Different Seeds at Same Centroid Point. Eps = 0.0003125")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Reference convergence line
    ref_x = np.array([min(seed1234_n), max(seed1234_n)])
    ref_y = seed1234_l1[0] * (ref_x / seed1234_l1[0])**(-0.5)
    ax1.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

    ax1.legend(title="Seeds", fontsize=9, title_fontsize=10)
    plt.tight_layout()
    plt.show()

def show_final_L1_error_bar_vs_epsilon():
    epsilon_to_final_row = {}

    # Find the last row for each epsilon group
    for i in range(len(data)):
        eps = epsilons[i]
        epsilon_to_final_row[eps] = i

    # Reverse-sorted epsilon values for visual clarity (big to small)
    sorted_epsilons = sorted(epsilon_to_final_row.keys(), reverse=True)
    final_l1_errors = [l1_errors[epsilon_to_final_row[eps]] for eps in sorted_epsilons]
    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    bar_positions = np.arange(len(sorted_epsilons))

    # Plot
    fig, ax = plt.subplots(figsize=(3, 4))
    ax.bar(bar_positions, final_l1_errors, color='steelblue')

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel("Epsilon Shell Value \n (log scale, decreasing)")
    ax.set_ylabel("Final L1 Error")
    ax.set_title("Final L1 Error vs Epsilon Shell Size")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()

def show_final_relative_error_bar_vs_epsilon():
    epsilon_to_final_row = {}

    # Find the last row for each epsilon group
    for i in range(len(data)):
        eps = epsilons[i]
        epsilon_to_final_row[eps] = i

    # Reverse-sorted epsilon values (largest to smallest)
    sorted_epsilons = sorted(epsilon_to_final_row.keys(), reverse=True)
    
    # Compute relative error (%) at final convergence point
    final_relative_errors = [
        (l1_errors[epsilon_to_final_row[eps]] / hardcoded_true_T_solution) * 100
        for eps in sorted_epsilons
    ]

    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    bar_positions = np.arange(len(sorted_epsilons))

    # Plot
    fig, ax = plt.subplots(figsize=(3, 4))
    ax.bar(bar_positions, final_relative_errors, color='crimson')

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel("Epsilon Shell Value \n (log scale, decreasing)")
    ax.set_ylabel("Final Relative Error (%)")
    ax.set_title("Final Relative Error vs Epsilon Shell Size")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()

## -- BOUNDARY --
def show_L1_error_vs_N_Walks_multiple_geometry_representations():
    # Load datasets for 3 interior points
    boundary_2020_data = np.loadtxt("../geometry_experiments/all_epsilon_convergence_2020.csv", delimiter=',', skiprows=1)
    boundary_4040_data = np.loadtxt("../geometry_experiments/all_epsilon_convergence_4040.csv", delimiter=',', skiprows=1)
    boundary_8080_data = np.loadtxt("../geometry_experiments/all_epsilon_convergence_8080.csv", delimiter=',', skiprows=1)

    # Unpack columns: assuming consistent structure [epsilon, n_walks, l1_error, time, steps, ...]
    boundary_2020_n = boundary_2020_data[:, 1]
    boundary_2020_l1 = boundary_2020_data[:, 2]

    boundary_4040_n = boundary_4040_data[:, 1]
    boundary_4040_l1 = boundary_4040_data[:, 2]
    
    boundary_8080_n = boundary_8080_data[:, 1]
    boundary_8080_l1 = boundary_8080_data[:, 2]

    # Plot
    fig, ax1 = plt.subplots(figsize=(7, 4))

    ax1.plot(boundary_2020_n, boundary_2020_l1, 'o-', linewidth=1, markersize=2, label="Boundary: 20 radius, 20 theta")
    ax1.plot(boundary_4040_n, boundary_4040_l1, 's-', linewidth=1, markersize=2, label="Boundary: 40 radius, 40 theta")
    ax1.plot(boundary_8080_n, boundary_8080_l1, 's-', linewidth=1, markersize=2, label="Boundary: 80 radius, 80 theta")

    # Log-Log Axes
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Number of Walks (N)")
    ax1.set_ylabel("Log L1 Error")
    ax1.set_title("Log-Log Convergence at Different Boundaries \n 2020, 4040, 8080 at Same Centroid Point. Eps = 0.0003125")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Reference convergence line
    ref_x = np.array([min(boundary_2020_n), max(boundary_2020_n)])
    ref_y = boundary_2020_l1[0] * (ref_x / boundary_2020_l1[0])**(-0.5)
    ax1.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

    ax1.legend(title="Seeds", fontsize=9, title_fontsize=10)
    plt.tight_layout()
    plt.show()

def show_relative_error_vs_N_Walks():
    fig, (ax1) = plt.subplots(nrows=2, figsize=(8, 10), sharex=True)
    
    for eps in unique_epsilons:
        mask = epsilons == eps
        
        x = n_walks[mask]
        l1_errors_ys = l1_errors[mask]
        y = 360.5597826843833
        
        
        # Plot on both subplots
        ax1.plot(x, y, marker='o', linestyle='-', linewidth=1, markersize=2, label=f"ε = {eps:.3g}")
    
    ## Log-Log Plot (Top)
    ax1.set_ylabel("Log L1 Error")
    ax1.set_xlabel("Number of Walks (N)")
    ax1.tick_params(labelbottom=True)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_title("Log-Log Convergence of WoSt at a Single Point")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    # Reference convergence line (O(1/√N))
    ref_x = np.array([min(n_walks), max(n_walks)])
    ref_y = l1_errors[0] * (ref_x / n_walks[0])**(-0.5)
    ax1.plot(ref_x, ref_y, 'k--', label="Reference: slope = -0.5", alpha=0.7)

    ## Display
    ax1.legend(title="Epsilon Values", fontsize=9, title_fontsize=10)
    # plt.title("Log-Log Plot of Convergence of Walk on Spheres\n at a Single Point for Different ε Values", fontsize=14)
    # plt.grid(True)
    plt.tight_layout()
    plt.show()

def show_full_and_zoomed_step_distribution():
    # Prepare step data grouped by epsilon
    epsilon_to_steps = {}
    for i in range(len(data)):
        eps = epsilons[i]
        step_count = steps[i]
        if eps not in epsilon_to_steps:
            epsilon_to_steps[eps] = []
        epsilon_to_steps[eps].append(step_count)

    # Sort and prepare data for plotting
    sorted_epsilons = sorted(epsilon_to_steps.keys(), reverse=True)
    step_distributions = [epsilon_to_steps[eps] for eps in sorted_epsilons]
    labels = [f"{eps:.0e}" for eps in sorted_epsilons]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    
    fig.suptitle("Steps per Walk vs Epsilon Shell Size", fontsize=14)

    # Full box plot
    ax1.boxplot(step_distributions, vert=True, patch_artist=True, labels=labels, showfliers=True)
    ax1.set_ylabel("Log Steps per Walk (Full Range)")
    ax1.set_yscale('log')
    ax1.set_title("Full Range (Log Scale)")
    ax1.grid(axis='y', linestyle='--', alpha=0.6)
    ax1.tick_params(labelbottom=True)

    # Zoomed box plot (focus on IQR / 0–20 region)
    ax2.boxplot(step_distributions, vert=True, patch_artist=True, labels=labels, showfliers=True)
    ax2.set_ylim(0, 20)
    ax2.set_ylabel("Steps per Walk")
    ax2.set_xlabel("Epsilon Shell Value (log scale, decreasing)")
    ax2.set_title("Zoomed In on IQR (0–20 Steps)")
    ax2.grid(axis='y', linestyle='--', alpha=0.6)

    ax2.set_title("Zoomed In on IQR (0–20 Steps)")
    
    # Caption-style explanation below the bottom axis
    # fig.text(0.5, 0.01, 
    #         "Smaller epsilon values increase precision but lead to longer walks.\n"
    #         "Top panel uses a log scale to show full distribution (including outliers).\n"
    #         "Bottom panel zooms into the interquartile range for clearer comparison.",
    #         ha='center', fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, bottom=0.2)
    plt.show()

def print_outlier_ratios():
    print(f"{'Epsilon':<10} {'Samples':<10} {'Q3':<8} {'>Q3 Count':<12} {'Outlier Ratio'}")
    print("-" * 55)

    sorted_epsilons = sorted(epsilon_to_steps.keys(), reverse=True)
    for eps in sorted_epsilons:
        step_counts = np.array(epsilon_to_steps[eps])
        q3 = np.percentile(step_counts, 75)
        count_above_q3 = np.sum(step_counts > q3)
        total = len(step_counts)
        ratio = count_above_q3 / total
        print(f"{eps:<10.0e} {total:<10} {q3:<8.1f} {count_above_q3:<12} {ratio:.3f}")

def show_outlier_ratio_per_epsilon():
    def print_outlier_ratios():
        print(f"{'Epsilon':<10} {'Samples':<10} {'Q3':<8} {'>Q3 Count':<12} {'Outlier Ratio'}")
        print("-" * 55)

        sorted_epsilons = sorted(epsilon_to_steps.keys(), reverse=True)
        for eps in sorted_epsilons:
            step_counts = np.array(epsilon_to_steps[eps])
            q3 = np.percentile(step_counts, 75)
            count_above_q3 = np.sum(step_counts > q3)
            total = len(step_counts)
            ratio = count_above_q3 / total
            print(f"{eps:<10.0e} {total:<10} {q3:<8.1f} {count_above_q3:<12} {ratio:.3f}")
    
    # Group steps by epsilon
    epsilon_to_steps = {}
    for i in range(len(data)):
        eps = epsilons[i]
        step_count = steps[i]
        if eps not in epsilon_to_steps:
            epsilon_to_steps[eps] = []
        epsilon_to_steps[eps].append(step_count)

    # Compute outlier ratios
    sorted_epsilons = sorted(epsilon_to_steps.keys(), reverse=True)
    outlier_ratios = []
    labels = [f"{eps:.0e}" for eps in sorted_epsilons]

    for eps in sorted_epsilons:
        step_counts = np.array(epsilon_to_steps[eps])
        q3 = np.percentile(step_counts, 75)
        num_outliers = np.sum(step_counts > q3)
        total = len(step_counts)
        ratio = num_outliers / total
        outlier_ratios.append(ratio)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    bar_positions = np.arange(len(sorted_epsilons))
    ax.bar(bar_positions, outlier_ratios, color='indianred')

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel("Epsilon Shell Value (log scale, decreasing)")
    ax.set_ylabel("Outlier Ratio (> Q3)")
    ax.set_title("Proportion of High-Step Outliers vs Epsilon Shell Size")

    fig.suptitle("Outlier Frequency Remains Stable — More Samples ≠ Higher Risk", fontsize=14)
    fig.text(
        0.5, 0.02,
        "Each bar shows the percentage of walks above the 75th percentile (Q3) for that epsilon.\n"
        "Although smaller ε values produce more outlier steps in raw numbers, the proportion of such walks remains roughly constant, \n indicating that outliers are an expected statistical consequence of \n larger sample sizes, not increased instability",
        ha='center', fontsize=10
    )
    
    print_outlier_ratios()
    
    

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()
    
def show_avg_steps_per_walk_vs_epsilon():
    # Group steps by epsilon
    epsilon_to_steps = {}
    for i in range(len(data)):
        eps = epsilons[i]
        step_count = steps[i]
        if eps not in epsilon_to_steps:
            epsilon_to_steps[eps] = []
        epsilon_to_steps[eps].append(step_count)

    # Sort epsilon groups (biggest to smallest)
    sorted_epsilons = sorted(epsilon_to_steps.keys(), reverse=True)
    avg_steps = [np.mean(epsilon_to_steps[eps]) for eps in sorted_epsilons]
    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    bar_positions = np.arange(len(sorted_epsilons))

    # Plot
    fig, ax = plt.subplots(figsize=(3, 4))
    ax.bar(bar_positions, avg_steps, color='darkgreen')

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel("Epsilon Shell Value \n (log scale, decreasing)")
    ax.set_ylabel("Average Steps per Walk")
    ax.set_title("Average Steps per Walk \n vs Epsilon Shell Size")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()
    
def show_avg_time_per_walk_vs_epsilon():
    epsilon_to_final_row = {}

    # Loop through the data and keep only the last row per epsilon
    for i in range(len(data)):
        eps = epsilons[i]
        epsilon_to_final_row[eps] = i  # overwrite, so we keep the last occurrence

    # Extract final cumulative time and n_walks for each epsilon group
    sorted_epsilons = sorted(epsilon_to_final_row.keys(), reverse=True)
    avg_time_per_walk = []
    for eps in sorted_epsilons:
        idx = epsilon_to_final_row[eps]
        final_time = cumulative_time[idx]
        total_walks = n_walks[idx]
        avg_time_per_walk.append(final_time / total_walks)

    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    
    # Convert seconds to milliseconds
    avg_time_per_walk_ms = [t * 1000 for t in avg_time_per_walk]

    # Plot
    fig, ax1 = plt.subplots(figsize=(3, 4))
    bar_positions = np.arange(len(sorted_epsilons))
    ax1.bar(bar_positions, avg_time_per_walk_ms, color='mediumslateblue')

    ax1.set_xticks(bar_positions)
    ax1.set_xticklabels(labels, rotation=45)
    ax1.set_xlabel("Epsilon Shell Value \n (log scale, decreasing)")
    ax1.set_ylabel("Average Time per Walk (ms)")
    ax1.set_title("Average Time per Walk \n vs Epsilon Shell Size")

    # fig.suptitle("Smaller ε Increases Computational Cost per Walk", fontsize=14)
    # fig.text(
    #     0.5, 0.02,
    #     "Final cumulative time divided by number of walks gives the average time per walk for each epsilon.\n"
    #     "As epsilon decreases, more steps are needed per walk, increasing runtime per estimate.",
    #     ha='center', fontsize=10
    # )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()

def show_avg_time_and_walks_vs_epsilon():
    epsilon_to_final_row = {}

    # Group by epsilon: keep last index per epsilon
    for i in range(len(data)):
        eps = epsilons[i]
        epsilon_to_final_row[eps] = i

    # Extract final values for each epsilon
    sorted_epsilons = sorted(epsilon_to_final_row.keys(), reverse=True)
    avg_time_per_walk_ms = []
    total_walks_per_epsilon = []

    for eps in sorted_epsilons:
        idx = epsilon_to_final_row[eps]
        final_time = cumulative_time[idx]
        total_walks = n_walks[idx]
        avg_time_per_walk_ms.append((final_time / total_walks) * 1000)  # convert to ms
        total_walks_per_epsilon.append(total_walks)

    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    bar_positions = np.arange(len(sorted_epsilons))

    # Plot
    fig, ax1 = plt.subplots(figsize=(10, 5))

    # Left Y-axis: Average time per walk (ms)
    ax1.bar(bar_positions, avg_time_per_walk_ms, color='mediumslateblue')
    ax1.set_ylabel("Avg Time per Walk (ms)", color='mediumslateblue')
    ax1.set_xlabel("Epsilon Shell Value (log scale, decreasing)")
    ax1.set_xticks(bar_positions)
    ax1.set_xticklabels(labels, rotation=45)
    ax1.tick_params(axis='y', labelcolor='mediumslateblue')

    # Right Y-axis: Total number of walks
    ax2 = ax1.twinx()
    ax2.plot(bar_positions, total_walks_per_epsilon, 'o--', color='darkorange', label='Total Walks')
    ax2.set_ylabel("Total Walks to Converge", color='darkorange')
    ax2.tick_params(axis='y', labelcolor='darkorange')

    # Titles and caption
    fig.suptitle("How Epsilon Affects Time per Walk and Total Walks", fontsize=14)
    ax1.set_title("Balancing Walk Cost vs Quantity as Epsilon Decreases")
    fig.text(
        0.5, 0.02,
        "Smaller epsilon shells require more steps per walk (↑ time/walk) and more walks for convergence.\n"
        "This chart compares both to understand which dominates total compute cost.",
        ha='center', fontsize=10
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()

def show_total_walks_vs_epsilon():
    epsilon_to_final_row = {}

    # Find the last row for each epsilon group
    for i in range(len(data)):
        eps = epsilons[i]
        epsilon_to_final_row[eps] = i

    # Extract total walks at convergence
    sorted_epsilons = sorted(epsilon_to_final_row.keys(), reverse=True)
    total_walks = [n_walks[epsilon_to_final_row[eps]] for eps in sorted_epsilons]
    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    bar_positions = np.arange(len(sorted_epsilons))

    # Plot
    fig, ax = plt.subplots(figsize=(3, 4))
    ax.bar(bar_positions, total_walks, color='darkorange')

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel("Epsilon Shell Value \n (log scale, decreasing)")
    ax.set_ylabel("Number of Walks to Converge")
    ax.set_title("N Walks Needed to Converge \n vs Epsilon Shell Size")

    # fig.suptitle("More Walks Are Needed as Epsilon Decreases", fontsize=14)
    # fig.text(
    #     0.5, 0.02,
    #     "Each bar shows how many Monte Carlo walks were performed before convergence was reached.\n"
    #     "Smaller epsilon values increase required precision, which leads to a higher number of total walks.",
    #     ha='center', fontsize=10
    # )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()

def show_cumulative_time_vs_epsilon():
    epsilon_to_final_row = {}

    # Find the last row for each epsilon group
    for i in range(len(data)):
        eps = epsilons[i]
        epsilon_to_final_row[eps] = i

    # Extract final cumulative execution time per epsilon
    sorted_epsilons = sorted(epsilon_to_final_row.keys(), reverse=True)
    cumulative_times = [cumulative_time[epsilon_to_final_row[eps]] for eps in sorted_epsilons]
    labels = [f"{eps:.0e}" for eps in sorted_epsilons]
    bar_positions = np.arange(len(sorted_epsilons))

    # Convert to milliseconds for better readability
    cumulative_times_ms = [t for t in cumulative_times]

    # Plot
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.bar(bar_positions, cumulative_times_ms, color='steelblue')

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel("Epsilon Shell Value (log scale, decreasing)")
    ax.set_ylabel("Time to Convergence (s)")
    ax.set_title("Time to Convergence vs Epsilon Shell Size")

    # fig.suptitle("Total Runtime Increases with Smaller Epsilon", fontsize=14)
    # fig.text(
    #     0.5, 0.02,
    #     "Final cumulative time per epsilon group shows total compute cost for convergence.\n"
    #     "As epsilon decreases, both walk length and number of walks increase total runtime.",
    #     ha='center', fontsize=10
    # )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plt.show()

if __name__ == "__main__":
    # show_L1_error_vs_time()
    # show_percent_error_vs_time()
    # show_l1_and_percent_error_vs_time()
    
    # show_L1_error_vs_N_Walks()
    # show_L1_error_vs_N_Walks_single()
    # show_L1_error_vs_N_Walks_multiple_points()
    # show_final_L1_error_bar_vs_epsilon()
    # show_final_relative_error_bar_vs_epsilon()
    
    # show_L1_error_vs_N_Walks_multiple_seeds()
    # show_L1_error_vs_N_Walks_multiple_geometry_representations()
    
    # show_time_vs_epsilon() ASS
    
    #  STEPS!!
    # show_full_and_zoomed_step_distribution()
    # show_outlier_ratio_per_epsilon()
    show_avg_steps_per_walk_vs_epsilon()
    
    # show_avg_time_and_walks_vs_epsilon() SUCKS
    
    # show_avg_time_per_walk_vs_epsilon() 
    # show_total_walks_vs_epsilon()
    # show_cumulative_time_vs_epsilon()
    