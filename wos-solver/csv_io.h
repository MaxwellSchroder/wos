#pragma once

#include "types.h"
#include <vector>

void readCSVandAppendSegmentsAdvancedWithType(const std::string& filename, std::vector<Polyline>& boundaryDirichlet, std::vector<Polyline>& boundaryNeumann, std::unordered_map<Vec2D,double>& boundaryTemperatureMap);
void readCSVandAppendSegments(const std::string& filename, std::vector<std::array<Vec2D,2>>& scene, PointCloud& cloud);
void writeResultsToCSV(const std::vector<ExperimentResult>& results, const std::string& filename);
void readInteriorPointsT(const std::string& filename, std::vector<std::tuple<Vec2D, double>>& points);