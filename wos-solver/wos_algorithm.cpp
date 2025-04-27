#include "wos_algorithm.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <fstream>

using namespace std;

float dot(Vec2D u, Vec2D v) { return real(conj(u)*v); }
float length( Vec2D u ) { return sqrt( norm(u) ); }
inline float getX(const Vec2D& v) { return std::real(v); }
inline float getY(const Vec2D& v) { return std::imag(v); }

// Create global_tree
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
float treeBasedTemperatureQuery( Vec2D x, const PointCloud& cloud ) {
    if (!global_tree) {
       cerr << "No global tree. Error out" << endl;
       return -1;
    }
 
    float query_arr[2] = {real(x), imag(x)};
 
    // re-initialise result every time
    resultSet.init(&nearest_index, &out_dist_sqr);
 
    // perform query for nearest neighbour
    global_tree->findNeighbors(resultSet, query_arr, nanoflann::SearchParameters(10));
 
    // cerr << "Nearest point is = ( " << global_cloud.pts[nearest_index].x << ", " << global_cloud.pts[nearest_index].y << ")" << endl;
    // cerr << "Temperature there  = " << global_cloud.pts[nearest_index].temperature << " Kelvin" << endl;
 
    return cloud.pts[nearest_index].temperature;
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

// function to set up all required point cloud function calls
void setupKDTree(const PointCloud& cloud) {
    // Determine the parameter maxLeaf
    size_t maxLeaf = std::max<size_t>(10,cloud.pts.size() / 10);
 
    // Initialise the global_tree
    global_tree = std::make_unique<KDTree>(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeaf));
 
    // build the tree
    global_tree->buildIndex();
}