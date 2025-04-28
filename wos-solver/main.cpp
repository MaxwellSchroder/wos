#include "types.h"
#include "wos_algorithm.h"
#include "experiment_runner.h"
#include "csv_io.h"
#include <iostream>

using namespace std;

int main() {
    // Read in the combined_coordinates, and generate the scene vector<Segment>
    vector<Segment> scene;

    PointCloud cloud;
    readCSVandAppendSegments("boundary_representation.csv", scene, cloud);
    setupKDTree(cloud);

    std::vector<std::tuple<Vec2D, float>> interior_points_T;
    readInteriorPointsT("interior_T_solution.csv", interior_points_T);

    // Testing a single point for convergence
    if (!interior_points_T.empty()) {
        int n_thetas = 20;
        int n_radii = 20;
        int flat_index = (n_thetas / 2) * n_radii + (n_radii / 2);
        std::cout << "Calculated flat index for middle point: " << flat_index << std::endl;

        auto [test_point, T_true] = interior_points_T[flat_index]; // test_point::(x,y), t_true::Int
        
        if (!insideDomain(test_point, scene)) {
            std::cerr << "WARNING: Selected test point is NOT inside the domain!\n";
        } else {
            std::cerr << "Point selected is inside of the domain!\n";
        }

        std::vector<float> epsilons = {0.05f, 0.02f, 0.01f, 0.005f, 0.002f, 0.001f, 0.0005f};
        const int nWalkLowerLimit = 1;
        const int nWalkUpperLimit = static_cast<int>(std::pow(2, 15));
        const int nWalkIncrement = 1;
        std::vector<ExperimentResult> results;

        for (float eps : epsilons) {
            // seed random for reproduceable results
            srand(1234); // srand( time(NULL) );

            std::cerr << "Solving for eps = " << eps << " ...";
            testSinglePointConvergence(test_point, scene, T_true, nWalkLowerLimit, nWalkUpperLimit, nWalkIncrement, eps, results, cloud);
        }

        writeResultsToCSV(results, "all_epsilon_convergence.csv");
    } else {
        std::cerr << "Failed to read interior points with temperature!" << std::endl;
    }
    
    std::cerr << "Finished!" << std::endl;
    return 0;
}