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

    return 0;
}