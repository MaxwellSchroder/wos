#pragma once

#include "types.h"
#include <functional>

float dot(Vec2D u, Vec2D v);
float length(Vec2D u);
inline float getX(const Vec2D& v);
inline float getY(const Vec2D& v);

float singleWalkEstimate(Vec2D x0, const std::vector<Segment>& segments, std::function<float(Vec2D)> g, float eps);
float treeBasedTemperatureQuery(Vec2D x, const PointCloud& cloud);
double signedAngle(Vec2D x, const std::vector<Segment>& P);
bool insideDomain(Vec2D x, const std::vector<Segment>& boundaryDirichlet);
void setupKDTree(const PointCloud& cloud);