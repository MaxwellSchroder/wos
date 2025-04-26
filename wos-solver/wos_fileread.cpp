// (Slow) implementation of Muller's 1956 Walk on Spheres algorithm
// Corresponds to the naïve estimator given in Equation 5 of
// Sawhney & Crane, Monte Carlo Geometry Processing (2020).
// NOTE: this code makes a few shortcuts for the sake of code brevity; may
// be more suitable for tutorials than for production code/evaluation.
// To compile: g++ -std=c++17 -O3 -pedantic -Wall -I./include wos_fileread.cpp -o wos
#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include "nanoflann.hpp"


using namespace std;
namespace std {
   template<>
   struct hash<std::complex<float>> {
       std::size_t operator()(const std::complex<float>& v) const {
           auto h1 = std::hash<float>{}(v.real());
           auto h2 = std::hash<float>{}(v.imag());
           return h1 ^ (h2 << 1);  // Combine the hashes
       }
   };

   template<>
   struct equal_to<std::complex<float>> {
       bool operator()(const std::complex<float>& lhs, const std::complex<float>& rhs) const {
           return std::abs(lhs.real() - rhs.real()) < 1e-6 &&
                  std::abs(lhs.imag() - rhs.imag()) < 1e-6;
       }
   };
}

// use std::complex to implement 2D vectors
using Vec2D = complex<float>;
float dot(Vec2D u, Vec2D v) { return real(conj(u)*v); }
float length( Vec2D u ) { return sqrt( norm(u) ); }
inline float getX(const Vec2D& v) { return std::real(v); }
inline float getY(const Vec2D& v) { return std::imag(v); }

// a segment is just a pair of points
using Segment = array<Vec2D,2>;

struct PointCloud {
   struct Point {
      float x, y, temperature;
   };

   vector<Point> pts;

   //must return the number of data points
   inline size_t kdtree_get_point_count() const { return pts.size(); }

   // returns the dim'th coordinate of the point
   inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
      if (dim == 0) return pts[idx].x;
      else return pts[idx].y;
   }

   template <class BBOX>
   bool kdtree_get_bbox(BBOX&) const {return false ;}
};

// Create global_cloud
PointCloud global_cloud;

// Declare a global pointer for the KDTree
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 2> KDTree;
std::unique_ptr<KDTree> global_tree;

// Initialise your query memory buffers first
size_t nearest_index;
float out_dist_sqr;
nanoflann::KNNResultSet<float> resultSet(1);


// returns the point on segment s closest to x
Vec2D closestPoint( Vec2D x, Segment s ) {
   Vec2D u = s[1]-s[0];
   float t = clamp(dot(x-s[0],u)/dot(u,u),0.f,1.f);
   return (1-t)*s[0] + t*s[1];
}

// returns a random value in the range [rMin,rMax]
float random( float rMin, float rMax ) {
   const float rRandMax = 1./(float)RAND_MAX;
   float u = rRandMax*(float)rand();
   return u*(rMax-rMin) + rMin;
}

// solves a Laplace equation Δu = 0 at x0, where the boundary is given
// by a collection of segments, and the boundary conditions are given
// by a function g that can be evaluated at any point in space
float solve( Vec2D x0, vector<Segment> segments, function<float(Vec2D)> g, int nWalks, float eps) {
   const int maxSteps = 128; // maximum walk length

   float sum = 0.;
   for( int i = 0; i < nWalks; i++ ) {
      Vec2D x = x0;
      float R;
      int steps = 0;
      do {
         R = numeric_limits<float>::max();
         for( auto s : segments ) {
            Vec2D p = closestPoint( x, s );
            R = min( R, length(x-p) );
         }
         float theta = random( 0., 2.*M_PI );
         x = x + Vec2D( R*cos(theta), R*sin(theta) );
         steps++;
      }
      while( R > eps && steps < maxSteps );

      sum += g(x);
   }
   return sum/nWalks; // Monte Carlo estimate
}

// This function performs a single walk and returns the estimate
float singleWalkEstimate(Vec2D x0, const std::vector<Segment>& segments, function<float(Vec2D)> g, float eps) {
   const int maxSteps = 128;
   Vec2D x = x0;
   float R;
   int steps = 0;
   do {
       R = numeric_limits<float>::max();
       for (auto& s : segments) {
           Vec2D p = closestPoint(x, s);
           R = min(R, length(x - p));
       }
       float theta = random(0., 2. * M_PI);
       x = x + Vec2D(R * cos(theta), R * sin(theta));
       steps++;
   } while (R > eps && steps < maxSteps);

   return g(x);
}

// search through the tree, and query the closest temperature value given
float treeBasedTemperatureQuery( Vec2D x ) {
   if (!global_tree) {
      cerr << "No global tree. Error out" << endl;
      return -1;
   }

   // work without vec2D first
   // Vec2D queryPoint = {0.2f, 0.1f};
   float query_arr[2] = {real(x), imag(x)};

   // re-initialise result every time
   resultSet.init(&nearest_index, &out_dist_sqr);

   // perform query for nearest neighbour
   global_tree->findNeighbors(resultSet, query_arr, nanoflann::SearchParameters(10));

   // cerr << "Nearest point is = ( " << global_cloud.pts[nearest_index].x << ", " << global_cloud.pts[nearest_index].y << ")" << endl;
   // cerr << "Temperature there  = " << global_cloud.pts[nearest_index].temperature << " Kelvin" << endl;

   return global_cloud.pts[nearest_index].temperature;
}

// these routines are not used by WoSt itself, but are rather used to check
// whether a given evaluation point is actually inside the domain
double signedAngle( Vec2D x, const vector<Segment>& P )
{
   double Theta = 0.;
   for( int i = 0; i < P.size(); i++ )
      // there is always two, so no need to run a secon dloop
      Theta += arg((P[i][1]-x)/(P[i][0]-x));
   return Theta;
}

// Returns true if the point x is contained in the region bounded by the Dirichlet
// and Neumann curves.  We assume these curves form a collection of closed polygons,
// and are given in a consistent counter-clockwise winding order.
bool insideDomain( Vec2D x,
                   const vector<Segment>& boundaryDirichlet)
{
   double Theta = signedAngle( x, boundaryDirichlet );
   const double delta = 1e-4; // numerical tolerance
   return abs(Theta-2.*M_PI) < delta; // boundary winds around x exactly once
}

// This function reads in the CSV file given from filenow, and takesn in the scene array and the pointcloud. It will add all the segments to the scene, and generate the point cloud
void readCSVandAppendSegments(const string& filename, vector<Segment>& scene) {
   ifstream file(filename);  // Open the file for reading. If it fails, error
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

   vector<Vec2D> points;
   string lineFromFile;
   
   while (std::getline(file, lineFromFile)) {
      stringstream ss(lineFromFile);  // Create a stringstream from the line
      string x_str, y_str, temperature_str;

      // Parse x and y values (comma-separated)
      if (getline(ss, x_str, ',') && getline(ss, y_str, ',') && getline(ss, temperature_str)) {
         float x = std::stof(x_str);  // Convert x string to float
         float y = std::stof(y_str);  // Convert y string to float
         float temperature = std::stof(temperature_str);

         // Create a Vec2D point from the parsed x and y values
         Vec2D point = Vec2D(x,y);

         // Add the point to the points array
         points.push_back(point);

         // ALSO: Add the x, y, and a Temperature value into the PointCloud
         global_cloud.pts.push_back({x, y, temperature});
      }
   }

   cout << "Number of points: " << points.size() << std::endl;

   // Reverse the points to get a +ve winding order
   std::reverse(points.begin(), points.end());

   // Now, create segments from consecutive points (A-B, B-C, etc.)
   for (size_t i = 1; i < points.size(); ++i) {
      Segment s = {{points[i - 1], points[i]}};
      scene.emplace_back(s);  // Create segments and append to scene
   }

   // Always add a closing segment to ensure the domain is closed
   if (!points.empty()) {
      Segment closing = {{points.back(), points.front()}};
      scene.emplace_back(closing);
      std::cerr << "Scene closed by adding final segment.\n";
   }
   std::cerr << "Scene now has " << scene.size() << " segments.\n";

   file.close();
}

// function to set up all required point cloud function calls
void setupKDTree() {
   // Determine the parameter maxLeaf
   size_t maxLeaf = std::max<size_t>(10,global_cloud.pts.size() / 10);

   // Initialise the global_tree
   global_tree = std::make_unique<KDTree>(2, global_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeaf));

   // build the tree
   global_tree->buildIndex();
}

void printScene(const std::vector<Segment>& scene3) {
   cerr << "Printing scene" << endl;
   for (const auto& segment : scene3) {
      cerr << "getting into loop" << endl;
      std::cerr << "Segment from (" << real(segment[0]) << ", " << imag(segment[0]) << ") "
               << "to (" << real(segment[1]) << ", " << imag(segment[1]) << ")\n";
   }
   cerr << "Done!" << endl;

}

void readInteriorPointsT(const std::string& filename, std::vector<std::tuple<Vec2D, float>>& points) {
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

           float x = std::stof(x_str);
           float y = std::stof(y_str);
           float t = std::stof(t_str);
           points.emplace_back(Vec2D(x, y), t);
       }
   }
}

std::unordered_map<Vec2D, double> readBoundaryTemperatureMap(const std::string& filename) {
   std::unordered_map<Vec2D, double> boundaryMap;

   std::ifstream file(filename);
   if (!file.is_open()) {
       std::cerr << "Error: could not open boundary CSV: " << filename << std::endl;
       return boundaryMap;
   }

   std::string line;
   while (std::getline(file, line)) {
       std::stringstream ss(line);
       std::string x_str, y_str, t_str;
       if (std::getline(ss, x_str, ',') &&
           std::getline(ss, y_str, ',') &&
           std::getline(ss, t_str, ',')) {

           float x = std::stof(x_str);
           float y = std::stof(y_str);
           double t = std::stod(t_str);

           Vec2D key(x, y);
           boundaryMap[key] = t;
       }
   }

   return boundaryMap;
}

// Outer function to cumulative grow estimates at a single point and record them to a CSV file
void testSinglePointConvergence(
   Vec2D x0,
   const std::vector<Segment>& scene,
   float T_true,
   int minWalks,
   int maxWalks,
   int walkCheckpointIncrement,
   float eps,
   const std::string& outputFile = "walk_vs_error.csv"
) {
   std::ofstream out(outputFile);
   if (!out.is_open()) {
       std::cerr << "Error: Could not open output file: " << outputFile << std::endl;
       return;
   }

   float running_sum = 0.0;
   int walks_completed = 0;

   for (int walk_counter = 1; walk_counter <= maxWalks; ++walk_counter) {
       float walk_result = singleWalkEstimate(x0, scene, treeBasedTemperatureQuery, eps);
       running_sum += walk_result;
       walks_completed++;

       // Output the error at checkpoint intervals
       if (walk_counter >= minWalks && (walk_counter - minWalks) % walkCheckpointIncrement == 0) {
           float T_estimate = running_sum / walks_completed;
           float l1_error = std::abs(T_estimate - T_true);

           out << walk_counter << "," << l1_error << "\n";
           std::cerr << "[Walks = " << walk_counter << "] T_est = " << T_estimate
                     << ", T_true = " << T_true << ", L1 error = " << l1_error << "\n";
       }
   }

   std::cerr << "Finished writing cumulative error data to " << outputFile << "\n";
}

void runInteriorEstimation(const std::vector<Vec2D>& interior_points,
                           const std::vector<Segment>& scene,
                           const std::unordered_map<Vec2D, double>& boundaryMap,
                           const std::string& outputFile) {
    using Entry = std::tuple<float, float, double>;
    std::vector<Entry> results;

    for (const auto& x0 : interior_points) {
        std::ostringstream key;
        key << std::fixed << std::setprecision(6) << getX(x0) << "," << getY(x0);

        double u;
        auto it = boundaryMap.find(x0);// check if point is already on the boundary
         if (it != boundaryMap.end()) {
            u = it->second;
         } else {
            u = solve(x0, scene, treeBasedTemperatureQuery, 128, 0.01);
         }

        results.emplace_back(getX(x0), getY(x0), u);
    }

    // Write results to CSV
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        std::cerr << "Error: could not open output file: " << outputFile << std::endl;
        return;
    }

    for (const auto& [x, y, t] : results) {
        out << x << "," << y << "," << t << "\n";
    }
}

int main( int argc, char** argv ) {
   // seed random for reproduceable results
   srand(1234);
   // srand( time(NULL) );

   // Read in the combined_coordinates, and generate the scene vector<Segment>
   vector<Segment> scene;
   readCSVandAppendSegments("boundary_representation.csv",scene);
   setupKDTree();
   auto boundaryMap = readBoundaryTemperatureMap("boundary_representation.csv"); // Data structure that stores "x,y" -> temperature values for O(1) lookup

   std::vector<std::tuple<Vec2D, float>> interior_points_T;
   readInteriorPointsT("interior_T_solution.csv", interior_points_T);

   // ofstream out( "out.csv" );

   // which technique will be used to solve the estimation
   // runInteriorEstimation(interior_points, scene, boundaryMap, "estimated_solution.csv"); // INTERIOR IS FOR EACH INTERIOR POINT

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

      const float eps = 0.01;
      const int nWalkLowerLimit = 1;
      const int nWalkUpperLimit = static_cast<int>(std::pow(2, 14));
      const int nWalkIncrement = 1;

      testSinglePointConvergence(test_point, scene, T_true, nWalkLowerLimit, nWalkUpperLimit, nWalkIncrement, eps, "error_plot_x0.csv");
   } else {
      std::cerr << "Failed to read interior points with temperature!" << std::endl;
   }

   std::cerr << "Finished!" << std::endl;
   return 0;
}
