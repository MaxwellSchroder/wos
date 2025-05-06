#include "csv_io.h"

#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

void readCSVandAppendSegmentsAdvancedWithType(
    const string& filename,
    vector<Polyline>& boundaryDirichlet,
    vector<Polyline>& boundaryNeumann,
    std::unordered_map<Vec2D,double>& boundaryTemperatureMap
) {
    ifstream file(filename);  // Open the file for reading. If it fails, error
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    vector<Vec2D> points;
    string lineFromFile;
    
    // Temporary polyline we're building
    Polyline currentPolyline;
    string currentType = "N";  // "Dirichlet" or "Neumann"
    
    while (std::getline(file, lineFromFile)) {
        stringstream ss(lineFromFile);  // Create a stringstream from the line
        string type_str, x_str, y_str, temperature_str;

        // Parse x and y values (comma-separated)
        if (getline(ss, type_str, ',') && getline(ss, x_str, ',') && getline(ss, y_str, ',') && getline(ss, temperature_str)) {
            double x = std::stof(x_str);  // Convert x string to double
            double y = std::stof(y_str);  // Convert y string to double
            double temperature = std::stof(temperature_str);

            // Create a Vec2D point from the parsed x and y values
            Vec2D point = Vec2D(x,y);

            // If neumann, add to Neumann. If Dirichlet, add to Dirichlet.
            if (currentType != type_str) {
                // Push old polyline onto its respective type
                if (currentType == "D") {
                    boundaryDirichlet.push_back({{currentPolyline}});
                } else {
                    boundaryNeumann.push_back({{currentPolyline}});
                }

                // Start of new Polyline
                currentPolyline.clear(); 
                currentType = type_str;
            }
            // Push this point onto the current Polyline
            currentPolyline.push_back(point);

            // If it's a Dirichlet, we should also add the (x,y) -> Temp hash mapping
            if (currentType == "D") {
                boundaryTemperatureMap[point] = temperature;
            }
        }
    }

    // Final push of current Polyline because it wouldn't have registered any changes
    if (currentType == "D") {
        boundaryDirichlet.push_back({{currentPolyline}});
    } else {
        boundaryNeumann.push_back({{currentPolyline}});
    }

    return;
}

// // This function reads in the CSV file given from filenow, and takesn in the scene array and the pointcloud. It will add all the segments to the scene, and generate the point cloud
// void readCSVandAppendSegments(const string& filename, vector<Segment>& scene, PointCloud& cloud) {
//     ifstream file(filename);  // Open the file for reading. If it fails, error
//      if (!file.is_open()) {
//          std::cerr << "Error opening file: " << filename << std::endl;
//          return;
//      }
 
//     vector<Vec2D> points;
//     string lineFromFile;
    
//     while (std::getline(file, lineFromFile)) {
//        stringstream ss(lineFromFile);  // Create a stringstream from the line
//        string x_str, y_str, temperature_str;
 
//        // Parse x and y values (comma-separated)
//        if (getline(ss, x_str, ',') && getline(ss, y_str, ',') && getline(ss, temperature_str)) {
//           double x = std::stof(x_str);  // Convert x string to double
//           double y = std::stof(y_str);  // Convert y string to double
//           double temperature = std::stof(temperature_str);
 
//           // Create a Vec2D point from the parsed x and y values
//           Vec2D point = Vec2D(x,y);
 
//           // Add the point to the points array
//           points.push_back(point);
 
//           // ALSO: Add the x, y, and a Temperature value into the PointCloud
//           cloud.pts.push_back({x, y, temperature});
//        }
//     }
 
//     cout << "Number of points: " << points.size() << std::endl;
 
//     // Reverse the points to get a +ve winding order
//     std::reverse(points.begin(), points.end());
 
//     // Now, create segments from consecutive points (A-B, B-C, etc.)
//     for (size_t i = 1; i < points.size(); ++i) {
//        Segment s = {{points[i - 1], points[i]}};
//        scene.emplace_back(s);  // Create segments and append to scene
//     }
 
//     // Always add a closing segment to ensure the domain is closed
//     if (!points.empty()) {
//        Segment closing = {{points.back(), points.front()}};
//        scene.emplace_back(closing);
//        std::cerr << "Scene closed by adding final segment.\n";
//     }
//     std::cerr << "Scene now has " << scene.size() << " segments.\n";
 
//     file.close();
// }

// This function takes out ExperimentResults struct, and then outputs to a known format in CSV
void writeResultsToCSV(const std::vector<ExperimentResult>& results, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << "\n";
        return;
    }
    out << "epsilon,nWalks,l1_error,cumulative_time\n"; // header
    for (const auto& res : results) {
        out << res.epsilon << "," << res.nWalks << "," << res.l1_error << "," << res.cumulative_time << "\n";
    }
    out.close();
}
 
void readInteriorPointsT(const std::string& filename, std::vector<std::tuple<Vec2D, double>>& points) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open interior points file: " << filename << std::endl;
        return;
    }
 
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string x_str, y_str, t_str;
        if (std::getline(ss, x_str, ',') &&
            std::getline(ss, y_str, ',') &&
            std::getline(ss, t_str, ',')) {
 
            double x = std::stof(x_str);
            double y = std::stof(y_str);
            double t = std::stof(t_str);
            points.emplace_back(Vec2D(x, y), t);
        }
    }
}
 