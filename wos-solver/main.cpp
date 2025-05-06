#include "types.h"
// #include "wos_algorithm.h"
#include "wost_algorithm.h"
#include "csv_io.h"
#include <iostream>
#include <chrono>

using namespace std;

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
 
    auto start = std::chrono::high_resolution_clock::now();
    for (int walk_counter = 1; walk_counter <= maxWalks; ++walk_counter) {
        double walk_result = singleWalkStarEstimate(x0, boundaryDirichlet, boundaryNeumann, g, eps);
        running_sum += walk_result;
        walks_completed++;
 
        // Output the error at checkpoint intervals
        if (walk_counter >= minWalks && (walk_counter - minWalks) % walkCheckpointIncrement == 0) {
            double T_estimate = running_sum / walks_completed;
            double l1_error = std::abs(T_estimate - T_true);

            // End timer and record time it took for specific walk
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            results.push_back({eps, walk_counter, l1_error, elapsed.count()});
        }
    }
    
    std::cerr << "WoST: Finished writing cumulative Walks and L1 Error data to results for Eps = " << eps << "\n";
}

vector<Polyline> boundaryDirichlet;
vector<Polyline> boundaryNeumann;

std::unordered_map<Vec2D, double> boundaryTemperatureMap;

// --- FUNCTION G to use to find the temperature for now ---
auto g (Vec2D x, Vec2D p0, Vec2D p1) -> double {
    try {
        double temp_p0 = boundaryTemperatureMap.at(p0);
        double temp_p1 = boundaryTemperatureMap.at(p1);

        std::cout << "[g(x)]: "
                  << "p0 = " << p0 << ", "
                  << "p1 = " << p1 << ", "
                  << "x = (" << real(x) << "," << imag(x) << "), " << "Sum of 2 Segment Temps = " << ((temp_p0 + temp_p1) / 2) << "\n";

        // For now just simple average (or later interpolation)

        // Compute interpolation factor t along the segment [p0, p1]
        double segment_length = std::abs(p1 - p0);
        if (segment_length < 1e-12) {
            std::cerr << "[g(x)] WARNING: Zero-length segment encountered.\n";
            return temp_p0;  // fallback to a single temperature value if they are this close
        }

        double t = std::abs(x - p0) / segment_length;
        t = std::clamp(t, 0.0, 1.0);  // safety clamp

        return (1.0 - t) * temp_p1 + t * temp_p0;
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
        std::vector<double> epsilons = {0.01, 0.005, 0.00125, 0.0005, 0.00025, 0.000125};
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
        }

        writeResultsToCSV(results, "all_epsilon_convergence.csv");

    } else {
        std::cerr << "Failed to read interior points with temperature!" << std::endl;
    }
    
    std::cerr << "Finished!" << std::endl;
    return 0;
}