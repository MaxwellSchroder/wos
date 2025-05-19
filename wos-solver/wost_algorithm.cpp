#include "wost_algorithm.h"

#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <optional>

using namespace std;

const double infinity = numeric_limits<double>::infinity();

// returns a random value in the range [rMin,rMax]
double random( double rMin, double rMax ) {
    const double rRandMax = 1.0/(double)RAND_MAX;
    double u = rRandMax*(double)rand();
    return u*(rMax-rMin) + rMin;
}

// --- Geometry and vector utilities ---
// using Vec2D = std::complex<double>;
double length( Vec2D u ) { return sqrt( norm(u) ); }
double angleOf(Vec2D u) { return arg(u); }
Vec2D rotate90( Vec2D u ) { return Vec2D( -imag(u), real(u) ); }
double   dot(Vec2D u, Vec2D v) { return real(u)*real(v) + imag(u)*imag(v); }
double cross(Vec2D u, Vec2D v) { return real(u)*imag(v) - imag(u)*real(v); }

// returns the closexst point to x on a segment with endpoints a and b
Vec2D closestPoint( Vec2D x, Vec2D a, Vec2D b ) {
    Vec2D u = b-a;
    double t = clamp( dot(x-a,u)/dot(u,u), 0.0, 1.0 );
    return (1.0-t)*a + t*b;
}

// --- Polyline and boundary operations ---

// returns true if the point b on the polyline abc is a silhoutte relative to x
bool isSilhouette( Vec2D x, Vec2D a, Vec2D b, Vec2D c ) {
    return cross(b-a,x-a) * cross(c-b,x-b) < 0;
}

// returns the time t at which the ray x+tv intersects segment ab,
// or infinity if there is no intersection
double rayIntersection( Vec2D x, Vec2D v, Vec2D a, Vec2D b ) {
    Vec2D u = b - a;
    Vec2D w = x - a;
    double d = cross(v,u);
    double s = cross(v,w) / d;
    double t = cross(u,w) / d;
    if (t > 0. && 0. <= s && s <= 1.) {
       return t;
    }
    return infinity;
}

// returns distance from x to closest point on the given polylines P
double distancePolylines( Vec2D x, const std::vector<Polyline>& P, Vec2D& hit_p0, Vec2D& hit_p1) {
    double d = infinity; // minimum distance so far
    for( size_t i = 0; i < P.size(); i++ ) { // iterate over polylines
        for( size_t j = 0; j < P[i].size()-1; j++ ) { // iterate over segments
            Vec2D y = closestPoint( x, P[i][j], P[i][j+1] ); // distance to segment

            // Add logic to ensure p0 is getting overwritten for correct segments
            if (d > length(x-y)) {
                hit_p0 = P[i][j];    // <-- Capture segment endpoints
                hit_p1 = P[i][j+1];
            }

            d = min( d, length(x-y) ); // update minimum distance
        }
    }
    return d;
}

// returns distance from x to closest silhouette point on the given polylines P
double silhouetteDistancePolylines( Vec2D x, const std::vector<Polyline>& P ){
    double d = infinity; // minimum distance so far
    for( size_t i = 0; i < P.size(); i++ ) { // iterate over polylines
        for( size_t j = 1; j < P[i].size()-1; j++ ) { // iterate over segment pairs
            if( isSilhouette( x, P[i][j-1], P[i][j], P[i][j+1] )) {
                d = min( d, length(x-P[i][j]) ); // update minimum distance
            }
        }
    }
    return d;
}

// finds the first intersection y of the ray x+tv with the given polylines P,
// restricted to a ball of radius r around x.  The flag onBoundary indicates
// whether the first hit is on a boundary segment (rather than the sphere), and
// if so sets n to the normal at the hit point.
Vec2D intersectPolylines(Vec2D x, Vec2D v, double r,
                         const vector<Polyline>& P,
                         Vec2D& n, bool& onBoundary) {
    double tMin = r; // smallest hit time so far
    n = Vec2D{ 0.0, 0.0 }; // first hit normal
    onBoundary = false; // will be true only if the first hit is on a segment
    for (size_t i = 0; i < P.size(); i++) { // iterate over polylines
        for (size_t j = 0; j < P[i].size() - 1; j++) { // iterate over segments
            const double c = 1e-5; // ray offset (to avoid self-intersection)
            double t = rayIntersection( x + c*v, v, P[i][j], P[i][j+1] );
            if( t < tMin ) { // closest hit so far
                tMin = t;
                n = rotate90( P[i][j+1] - P[i][j] ); // get normal
                n /= length(n); // make normal unit length
                onBoundary = true;
            }
        }
    }
    return x + tMin * v; // first hit location
}

// --- Walk on Stars Solver ---
std::optional<double> singleWalkStarEstimate(
    Vec2D x0,
    const vector<Polyline>& boundaryDirichlet,
    const vector<Polyline>& boundaryNeumann,
    std::function<double(Vec2D x, Vec2D p0, Vec2D p1)> g,
    double eps
) {
    const double rMin = 0.0001;
    const int maxSteps = 65536; // typical for single walks

    Vec2D x = x0; // start walk at the evaluation point
    Vec2D n{ 0.0, 0.0 }; // assume x0 is an interior point, and has no normal
    bool onBoundary = false; // flag whether x is on the interior or boundary

    Vec2D hit_p0, hit_p1; // to catch the segment in which the polylines intersected

    int steps = 0;
    double r, dDirichlet, dSilhouette; // radii used to define star shaped region
    do { // loop until the walk hits the Dirichlet boundary
        // Compute star radius
        dDirichlet = distancePolylines(x, boundaryDirichlet, hit_p0, hit_p1);
        dSilhouette = silhouetteDistancePolylines(x, boundaryNeumann);
        r = max(rMin, min(dDirichlet, dSilhouette));

        // intersect a ray with the star-shaped region boundary
        double theta = random( -M_PI, M_PI );
        if (onBoundary) { // sample from a hemisphere around the normal
            theta = theta / 2.0 + angleOf(n);
        }
        Vec2D v{ cos(theta), sin(theta) }; // unit ray direction

        // Move to intersection point
        x = intersectPolylines(x, v, r, boundaryNeumann, n, onBoundary);

        steps++;
    }
    while (dDirichlet > eps && steps < maxSteps);

    if (steps >= maxSteps || std::isnan(x.real()) || std::isnan(x.imag())) {
        std::cerr << "[Walk Failed] Reached max steps or invalid point.\n";
        return std::nullopt;
    }

    return g(x, hit_p0, hit_p1);
}

// --- Domain Checking ---

// these routines are not used by WoSt itself, but are rather used to check
// whether a given evaluation point is actually inside the domain
double signedAngle( Vec2D x, const vector<Polyline>& P ) 
{
    double Theta = 0.;
    for( size_t i = 0; i < P.size(); i++ )
        for( size_t j = 0; j < P[i].size()-1; j++ )
            Theta += arg( (P[i][j+1]-x)/(P[i][j]-x) );
    return Theta;
}

// Returns true if the point x is contained in the region bounded by the Dirichlet
// and Neumann curves.  We assume these curves form a collection of closed polygons,
// and are given in a consistent counter-clockwise winding order.
bool insideDomain( Vec2D x,
    const vector<Polyline>& boundaryDirichlet,
    const vector<Polyline>& boundaryNeumann ) 
{
    double Theta = signedAngle( x, boundaryDirichlet ) +signedAngle( x, boundaryNeumann );
    const double delta = 1e-4; // numerical tolerance
    return abs(Theta-2.*M_PI) < delta; // boundary winds around x exactly once
}