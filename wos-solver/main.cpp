#include "types.h"
// #include "wos_algorithm.h"
#include "wost_algorithm.h"
#include "csv_io.h"
#include <iostream>
#include <chrono>
#include <tuple>

using namespace std;

std::pair<double, double> computeMeanAndStdDev(const std::deque<double>& values) {
    if (values.empty()) {
        return {0.0, 0.0};
    }

    double sum = 0.0;
    for (double v : values) {
        sum += v;
    }
    double mean = sum / values.size();

    double sq_sum = 0.0;
    for (double v : values) {
        sq_sum += (v - mean) * (v - mean);
    }

    double stddev = 0.0;
    if (values.size() > 1) {
        stddev = std::sqrt(sq_sum / (values.size() - 1));  // unbiased sample stddev
    }

    return {mean, stddev};
}

// --- testSinglePointUsingStars ---
void testSinglePointConvergenceWoSt(
    Vec2D x0,
    vector<Polyline> boundaryDirichlet, // absorbing part of the boundary
    vector<Polyline> boundaryNeumann, // reflecting part of boudnary
    double T_true,
    int minWalks,
    int maxWalks,
    int walkCheckpointIncrement,
    double eps,
    std::vector<ExperimentResult>& results,
    std::function<double(Vec2D x, Vec2D p0, Vec2D p1)> g
) {
    double running_sum = 0.0;
    int walks_completed = 0;
    int failed_walk_count = 0;
    int min_walks_required = static_cast<int>(std::pow(2, 10));

    std::deque<double> recent_l1_errors;
    int rse_buffer_size = 2000;
    const double rse_threshold = 0.01;
    std::cerr << "Ignoring maxwalks" << maxWalks << "\n";
    
    auto start = std::chrono::high_resolution_clock::now();

    int walk_counter = 1;
    while (true) {
        rse_buffer_size = std::max(2000, std::min(walk_counter / 10, 10000));

        // You have a minimum number of walks && RSE check if buffer is full
        if ((walk_counter > min_walks_required) && (recent_l1_errors.size() == static_cast<size_t>(rse_buffer_size)) && (walk_counter%rse_buffer_size == 0)) {
            auto [mean, stddev] = computeMeanAndStdDev(recent_l1_errors);
            double rse = std::abs(stddev / mean);

            std::cerr << "[RSE Check] ε = " << eps << ", RSE = " << rse << "\n";

            if (rse < rse_threshold) {
                std::cerr << "[Convergence] RSE threshold met. ε = " << eps << ", Walks = " << walk_counter << "\n";
                break;
            } else {
                std::cerr << "[Convergence NOT met] RSE threshold not met. ε = " << eps << ", Walks = " << walk_counter << "\n";
            }
        }
        
        double walk_result;
        int steps_in_walk;
        // Just in case walk hits step limit
        while (true) {
            auto maybe_result = singleWalkStarEstimate(x0, boundaryDirichlet, boundaryNeumann, g, eps);
            
            if (maybe_result.has_value()) {
                std::tie(walk_result, steps_in_walk) = maybe_result.value(); // or *maybe_result
                break;
                // use result and steps_in_walk
            } else {
                // handle failure case
                failed_walk_count++;
                std::cerr << "Retrying walk... (total failures so far: " << failed_walk_count << ")\n";
            }
        }

        running_sum += walk_result;
        walks_completed++;

        // Output the error at checkpoint intervals
        if (walk_counter >= minWalks && (walk_counter - minWalks) % walkCheckpointIncrement == 0) {

            double T_estimate = running_sum / walks_completed;

            double l1_error = std::abs(T_estimate - T_true);

            if (std::isnan(l1_error) || std::isinf(l1_error)) {
                l1_error = 0.0;  // treat as no error if both are extremely close
            }

            // End timer and record time it took for specific walk
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            results.push_back({eps, walk_counter, l1_error, elapsed.count(), steps_in_walk});

            // Update rolling buffer
            recent_l1_errors.push_back(l1_error);
            if (recent_l1_errors.size() > static_cast<size_t>(rse_buffer_size)) {
                recent_l1_errors.pop_front();
            }
        }
        // Manual increment of walk counter
        walk_counter++;
        
    }
    std::cerr << "WoST: Finished for Eps = " << eps << " | Failed walk count (number of exceeded maxWalks): " << failed_walk_count << " | Valid walks: " << walks_completed << "\n";
}

vector<Polyline> boundaryDirichlet;
vector<Polyline> boundaryNeumann;

std::unordered_map<Vec2D, double> boundaryTemperatureMap;

// --- FUNCTION G to use to find the temperature for now ---
auto g (Vec2D x, Vec2D p0, Vec2D p1) -> double {
    try {
        if (std::isnan(real(p0)) || std::isnan(imag(p0)) ||
            std::isnan(real(p1)) || std::isnan(imag(p1)) ||
            std::isnan(real(x)) || std::isnan(imag(x))) {
            std::cerr << "[g(x)] ERROR: NaN found in input coordinates.\n";
            std::cerr << "x = (" << real(x) << "," << imag(x) << "), "
                      << "p0 = (" << real(p0) << "," << imag(p0) << "), "
                      << "p1 = (" << real(p1) << "," << imag(p1) << ")\n";
            return 0.0;
        }

        double temp_p0 = boundaryTemperatureMap.at(p0);
        double temp_p1 = boundaryTemperatureMap.at(p1);

        if (std::isnan(temp_p0) || std::isnan(temp_p1)) {
            std::cerr << "[g(x)] ERROR: Temp lookup returned NaN.\n";
            std::cerr << "temp_p0 = " << temp_p0 << ", temp_p1 = " << temp_p1 << "\n";
            return 0.0;
        }

        // Compute interpolation factor t along the segment [p0, p1]
        double segment_length = std::abs(p1 - p0);
        if (segment_length < 1e-12) {
            std::cerr << "[g(x)] WARNING: Zero-length segment encountered.\n";
            return temp_p0;  // fallback to a single temperature value if they are this close
        }

        double t = std::abs(x - p0) / segment_length;
        t = std::clamp(t, 0.0, 1.0);  // safety clamp

        double interpolated_T = (1.0 - t) * temp_p1 + t * temp_p0;

        if (std::isnan(interpolated_T) || std::isinf(interpolated_T)) {
            std::cerr << "[g(x)] ERROR: Interpolated temperature is invalid.\n";
        }

        return interpolated_T;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "[g(x)] ERROR: Vertex not found in boundaryTemperatureMap.\n";
        std::cerr << "Attempted p0 = (" << real(p0) << "," << imag(p0) << "), "
                  << "p1 = (" << real(p1) << "," << imag(p1) << "), "
                  << "walk x = (" << real(x) << "," << imag(x) << ")\n";
        std::cerr << "Exception: " << e.what() << "\n";
        return 0.0; // Safe fallback
    }
};

double estimateInstrumentationOverheadPerWalk(int nMaxWalk = 10000) {
    using namespace std::chrono;

    std::vector<ExperimentResult> dummy_results;
    double dummy_eps = -1;
    int dummy_n = -1;
    double dummy_error = -1;

    auto start = high_resolution_clock::now();
    for (int i = 0; i < nMaxWalk; ++i) {
        auto t0 = high_resolution_clock::now();
        auto t1 = high_resolution_clock::now();
        duration<double> elapsed = t1 - t0;
        dummy_results.push_back({dummy_eps, dummy_n, dummy_error, elapsed.count(), 2});
    }
    auto end = high_resolution_clock::now();
    duration<double> total = end - start;

    return total.count() / nMaxWalk;
}

void print_boundaries() {
    // for simplicity, in this code we assume that the Dirichlet and Neumann
    // boundary polylines form a collection of closed polygons (possibly with holes),
    // and are given with consistent counter-clockwise orientation
    // vector<Polyline> boundaryDirichlet = {
    //     {{ Vec2D(0.2, 0.2), Vec2D(0.6, 0.0), Vec2D(1.0, 0.2) }},
    //     {{ Vec2D(1.0, 1.0), Vec2D(0.6, 0.8), Vec2D(0.2, 1.0) }}
    // };

    // vector<Polyline> boundaryNeumann = {
    //     {{ Vec2D(1.0, 0.2), Vec2D(0.8, 0.6), Vec2D(1.0, 1.0) }},
    //     {{ Vec2D(0.2, 1.0), Vec2D(0.0, 0.6), Vec2D(0.2, 0.2) }}
    // };
    for (const auto& [pos, temp] : boundaryTemperatureMap) {
        std::cout << "(" << pos.real() << ", " << pos.imag() << ") => " << temp << "\n";
    }

    std::cout << "boundaryDirichlet:" << std::endl;
    for (size_t i = 0; i < boundaryDirichlet.size(); ++i) {
        std::cout << "  Polyline " << i << ": { ";
        for (const auto& pt : boundaryDirichlet[i]) {
            std::cout << "(" << pt.real() << ", " << pt.imag() << ") ";
        }
        std::cout << "}" << std::endl;
    }

    std::cout << "boundaryNeumann:" << std::endl;
    for (size_t i = 0; i < boundaryNeumann.size(); ++i) {
        std::cout << "  Polyline " << i << ": { ";
        for (const auto& pt : boundaryNeumann[i]) {
            std::cout << "(" << pt.real() << ", " << pt.imag() << ") ";
        }
        std::cout << "}" << std::endl;
    }
}

int main() {
    // Attempt to read in boundary, both dirichlet and Neumann. Generate the Polylines.
    readCSVandAppendSegmentsAdvancedWithType("boundary_representation.csv", boundaryDirichlet, boundaryNeumann, boundaryTemperatureMap);

    // --- Read in boundary of dirichlet and Neumann
    std::vector<std::tuple<Vec2D, double>> interior_points_T;
    readInteriorPointsT("interior_T_solution.csv", interior_points_T);

    // Testing a single point for convergence
    if (!interior_points_T.empty()) {
        int n_thetas = 20;
        int n_radii = 20;
        int flat_index = (n_thetas / 2) * n_radii + (n_radii / 2);
        std::cout << "Calculated flat index for middle point: " << flat_index << std::endl;

        auto [test_point, T_true] = interior_points_T[flat_index]; // test_point::(x,y), t_true::Int
        std::cout << "test_point = " << test_point << " and T_true" << T_true;
        
        // std::vector<double> epsilons = {0.01, 0.005, 0.00125, 0.0005, 0.00025, 0.000125, 5e-05, 2.5e-05, 1.25e-05};
        std::vector<double> epsilons = {0.005, 0.00125, 0.0005, 0.00025, 0.000125};
        // std::vector<double> epsilons = {0.005};
        const int nWalkLowerLimit = 1;
        const int nWalkUpperLimit = static_cast<int>(std::pow(2, 15));
        const int nWalkIncrement = 1;
        std::vector<ExperimentResult> results;

        for (double eps : epsilons) {
            // seed random for reproduceable results
            srand(1234); // srand( time(NULL) );

            std::cerr << "Solving for eps = " << eps << " ...\n";


            testSinglePointConvergenceWoSt(test_point, boundaryDirichlet, boundaryNeumann, T_true, nWalkLowerLimit, nWalkUpperLimit, nWalkIncrement, eps, results, g);

        
            #ifdef ENABLE_INSTRUMENTATION
            double instrumentation_overhead = estimateInstrumentationOverheadPerWalk(nWalkUpperLimit);

             // Adjust cumulative_time for all result entries with current eps
            for (auto& r : results) {
                if (r.epsilon == eps) { // match current epsilon
                    double adjusted_time = r.cumulative_time - instrumentation_overhead*r.nWalks;
                    r.cumulative_time = std::max(0.0, adjusted_time); // clamp to non-negative
                }
            }
            #endif

            // Print final convergence time for this epsilon
            auto it = std::find_if(results.rbegin(), results.rend(), [eps](const ExperimentResult& r) {
                return r.epsilon == eps;
            });

            if (it != results.rend()) {
                std::cout << "[Convergence] ε = " << eps 
                        << " converged in " << it->cumulative_time << " seconds "
                        << "after " << it->nWalks << " walks.\n";
            } else {
                std::cout << "[Convergence] ε = " << eps << " produced no results.\n";
            }
        }

        writeResultsToCSV(results, "all_epsilon_convergence.csv");

    } else {
        std::cerr << "Failed to read interior points with temperature!" << std::endl;
    }
    
    std::cerr << "Finished!" << std::endl;
    return 0;
}