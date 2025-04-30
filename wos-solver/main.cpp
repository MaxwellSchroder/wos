#include "types.h"
// #include "wos_algorithm.h"
#include "wost_algorithm.h"
#include "csv_io.h"
#include <iostream>

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
 
    for (int walk_counter = 1; walk_counter <= maxWalks; ++walk_counter) {
        double walk_result = singleWalkStarEstimate(x0, boundaryDirichlet, boundaryNeumann, g, eps);
        running_sum += walk_result;
        walks_completed++;
 
        // Output the error at checkpoint intervals
        if (walk_counter >= minWalks && (walk_counter - minWalks) % walkCheckpointIncrement == 0) {
            double T_estimate = running_sum / walks_completed;
            double l1_error = std::abs(T_estimate - T_true);
            results.push_back({eps, walk_counter, l1_error});

            // std::cerr << "[Walks = " << walk_counter << "] T_est = " << T_estimate
            //           << ", T_true = " << T_true << ", L1 error = " << l1_error << "\n";
        }
    }
    
    std::cerr << "WoST: Finished writing cumulative Walks and L1 Error data to results for Eps = " << eps << "\n";
}


// for simplicity, in this code we assume that the Dirichlet and Neumann
// boundary polylines form a collection of closed polygons (possibly with holes),
// and are given with consistent counter-clockwise orientation
vector<Polyline> boundaryDirichlet = {
    {{ Vec2D(0.2, 0.2), Vec2D(0.6, 0.0), Vec2D(1.0, 0.2) }},
    {{ Vec2D(1.0, 1.0), Vec2D(0.6, 0.8), Vec2D(0.2, 1.0) }}
};

std::unordered_map<Vec2D, double> boundaryTemperatureMap;

vector<Polyline> boundaryNeumann = {
    {{ Vec2D(1.0, 0.2), Vec2D(0.8, 0.6), Vec2D(1.0, 1.0) }},
    {{ Vec2D(0.2, 1.0), Vec2D(0.0, 0.6), Vec2D(0.2, 0.2) }}
};

// --- FUNCTION G to use to find the temperature for now ---
auto g (Vec2D x, Vec2D p0, Vec2D p1) -> double {
    try {
        double temp_p0 = boundaryTemperatureMap.at(p0);
        double temp_p1 = boundaryTemperatureMap.at(p1);

        std::cout << "[g(x)] Successful lookup: "
                  << "p0 = " << p0 << ", "
                  << "p1 = " << p1 << ", "
                  << "x = (" << real(x) << "," << imag(x) << ")\n";

        // For now just simple average (or later interpolation)
        return (temp_p0 + temp_p1) / 2;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "[g(x)] ERROR: Vertex not found in boundaryTemperatureMap.\n";
        std::cerr << "Attempted p0 = (" << real(p0) << "," << imag(p0) << "), "
                  << "p1 = (" << real(p1) << "," << imag(p1) << "), "
                  << "walk x = (" << real(x) << "," << imag(x) << ")\n";
        std::cerr << "Exception: " << e.what() << "\n";
        return 0.0; // Safe fallback
    }

    // const double eps = 1e-2;

    // // Example: assign a constant temperature per polyline
    // // Here we make it 200 + 50*i Kelvin for polyline i
    // const double baseTemp = 200.0;
    // const double tempIncrement = 50.0;

    // for (size_t i = 0; i < boundaryDirichlet.size(); ++i) {
    //     const Polyline& poly = boundaryDirichlet[i];

    //     for (size_t j = 0; j + 1 < poly.size(); ++j) {
    //         Vec2D p0 = poly[j];
    //         Vec2D p1 = poly[j+1];

    //         double d0 = std::abs(x - p0);
    //         double d1 = std::abs(x - p1);
    //         double dSeg = std::abs(p1 - p0);

    //         if (std::abs((d0 + d1) - dSeg) < eps) {
    //             // Found the segment where x lies between p0 and p1
    //             std::cout << "Intersection between " << p0 << " and " << p1 << ". \n" << std::endl;
    //             return baseTemp + tempIncrement * i;  
    //             // 200K for polyline 0, 250K for polyline 1, etc.
    //         }
    //     }
    // }
};

int main() {
    boundaryTemperatureMap[Vec2D(0.2, 0.2)] = 200.0;
    boundaryTemperatureMap[Vec2D(0.6, 0.0)] = 200.0;
    boundaryTemperatureMap[Vec2D(1.0, 0.2)] = 200.0;

    boundaryTemperatureMap[Vec2D(1.0, 1.0)] = 450.0;
    boundaryTemperatureMap[Vec2D(0.6, 0.8)] = 450.0;
    boundaryTemperatureMap[Vec2D(0.2, 1.0)] = 450.0;

    // --- Read in boundary of dirichlet and Neumann
    std::vector<std::tuple<Vec2D, double>> interior_points_T;
    readInteriorPointsT("interior_T_solution.csv", interior_points_T);

    // Testing a single point for convergence
    if (!interior_points_T.empty()) {
        // int n_thetas = 20;
        // int n_radii = 20;
        // int flat_index = (n_thetas / 2) * n_radii + (n_radii / 2);
        // std::cout << "Calculated flat index for middle point: " << flat_index << std::endl;

        // auto [test_point, T_true] = interior_points_T[flat_index]; // test_point::(x,y), t_true::Int
        Vec2D test_point = Vec2D(0.6, 0.6);
        double T_true = 225;
        
        // std::vector<double> epsilons = {0.05f, 0.02f, 0.01f, 0.005f, 0.002f, 0.001f, 0.0005f};
        std::vector<double> epsilons = {0.05f};
        const int nWalkLowerLimit = 1;
        const int nWalkUpperLimit = static_cast<int>(std::pow(2, 10));
        // const int nWalkUpperLimit = 5;
        const int nWalkIncrement = 1;
        std::vector<ExperimentResult> results;

        for (double eps : epsilons) {
            // seed random for reproduceable results
            srand(1234); // srand( time(NULL) );

            std::cerr << "Solving for eps = " << eps << " ...";
            testSinglePointConvergenceWoSt(test_point, boundaryDirichlet, boundaryNeumann, T_true, nWalkLowerLimit, nWalkUpperLimit, nWalkIncrement, eps, results, g);
        }

        writeResultsToCSV(results, "all_epsilon_convergence.csv");
    } else {
        std::cerr << "Failed to read interior points with temperature!" << std::endl;
    }
    
    std::cerr << "Finished!" << std::endl;
    return 0;
}