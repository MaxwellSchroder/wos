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
#include <nanoflann.hpp>

using namespace std;

// use std::complex to implement 2D vectors
using Vec2D = complex<float>;
float dot(Vec2D u, Vec2D v) { return real(conj(u)*v); }
float length( Vec2D u ) { return sqrt( norm(u) ); }

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
float solve( Vec2D x0, vector<Segment> segments, function<float(Vec2D)> g ) {
   const float eps = 0.01; // stopping tolerance
   const int nWalks = 128; // number of Monte Carlo samples 4096
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

   // Now, create segments from consecutive points (A-B, B-C, etc.)
   for (size_t i = 1; i < points.size(); ++i) {
      Segment s = {{points[i - 1], points[i]}};
      scene.emplace_back(s);  // Create segments and append to scene
   }

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
   cerr << "DOne!" << endl;

}

int main( int argc, char** argv ) {
   // Read in the combined_coordinates, and generate the scene vector<Segment>
   vector<Segment> scene;

   readCSVandAppendSegments("combined_coordinates.csv",scene);
   setupKDTree();

   srand( time(NULL) );
   ofstream out( "out.csv" );

   int s = 256; // image size
   for( int j = 0; j < s; j++ )
   {
      cerr << "row "  << j <<  " of " << s << endl;
      for( int i = 0; i < s; i++ )
      {
         Vec2D x0( (float)i/(float)s, (float)j/(float)s );
         // double u = solve( x0, scene3, customTemperature );
         double u = 0.;
         // check that the point x0 in the image is actually inside of the domain.
         if (insideDomain(x0, scene)) {
            u = solve( x0, scene, treeBasedTemperatureQuery );
         }
         
         out << u;
         if( i < s-1 ) out << ",";
      }
      out << endl;
   }
   return 0;
   cerr << "Finished!" << endl;
}
