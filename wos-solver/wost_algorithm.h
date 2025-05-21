#pragma once

#include "types.h"
#include <vector>
#include <functional>
#include <complex>

// --- Geometry and vector utilities ---

double dot(Vec2D u, Vec2D v);
double length(Vec2D u);
double cross(Vec2D u, Vec2D v);
double angleOf(Vec2D u);
double random(double rMin, double rMax);
Vec2D rotate90(Vec2D u);

// --- Polyline and boundary operations ---
Vec2D closestPoint( Vec2D x, Vec2D a, Vec2D b );
bool isSilhouette(Vec2D x, Vec2D a, Vec2D b, Vec2D c);
double rayIntersection(Vec2D x, Vec2D v, Vec2D a, Vec2D b);
double distancePolylines(Vec2D x, const std::vector<Polyline>& P, Vec2D& hit_p0, Vec2D& hit_p1);
double silhouetteDistancePolylines(Vec2D x, const std::vector<Polyline>& P);
Vec2D intersectPolylines(Vec2D x, Vec2D v, double r, const std::vector<Polyline>& P, Vec2D& n, bool& onBoundary);

// --- Walk on Stars Solver ---

std::optional<std::tuple<double, float>> singleWalkStarEstimate(
    Vec2D x0,
    const std::vector<Polyline>& boundaryDirichlet,
    const std::vector<Polyline>& boundaryNeumann,
    std::function<double(Vec2D x, Vec2D p0, Vec2D p1)> g,
    double eps);

// --- Domain Checking ---

double signedAngle( Vec2D x, const std::vector<Polyline>& P );
bool insideDomain( Vec2D x,
    const std::vector<Polyline>& boundaryDirichlet,
    const std::vector<Polyline>& boundaryNeumann );