#pragma once

#include "types.h"
#include <vector>

void readCSVandAppendSegments(const std::string& filename, std::vector<std::array<Vec2D,2>>& scene);
void writeResultsToCSV(const std::vector<ExperimentResult>& results, const std::string& filename);
void readInteriorPointsT(const std::string& filename, std::vector<std::tuple<Vec2D, float>>& points);