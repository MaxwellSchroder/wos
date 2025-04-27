#include "experiment_runner.h"
#include <iostream>

using namespace std;

// Outer function to cumulative grow estimates at a single point and record them to a CSV file
void testSinglePointConvergence(
    Vec2D x0,
    const std::vector<Segment>& scene,
    float T_true,
    int minWalks,
    int maxWalks,
    int walkCheckpointIncrement,
    float eps,
    std::vector<ExperimentResult>& results,
    const PointCloud& cloud
 ) {
    auto g = [&cloud](Vec2D x) {
        return treeBasedTemperatureQuery(x, cloud);
    };

    float running_sum = 0.0;
    int walks_completed = 0;
 
    for (int walk_counter = 1; walk_counter <= maxWalks; ++walk_counter) {
        float walk_result = singleWalkEstimate(x0, scene, g, eps);
        running_sum += walk_result;
        walks_completed++;
 
        // Output the error at checkpoint intervals
        if (walk_counter >= minWalks && (walk_counter - minWalks) % walkCheckpointIncrement == 0) {
            float T_estimate = running_sum / walks_completed;
            float l1_error = std::abs(T_estimate - T_true);
            results.push_back({eps, walk_counter, l1_error});
 
          //   std::cerr << "[Walks = " << walk_counter << "] T_est = " << T_estimate
          //             << ", T_true = " << T_true << ", L1 error = " << l1_error << "\n";
        }
    }
    
    std::cerr << "Finished writing cumulative Walks and L1 Error data to results for Eps = " << eps << "\n";
}