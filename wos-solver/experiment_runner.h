#pragma once

#include "types.h"
#include "wos_algorithm.h"
#include <vector>

void testSinglePointConvergence(Vec2D x0, const std::vector<std::array<Vec2D,2>>& scene,
                                 float T_true, int minWalks, int maxWalks,
                                 int walkCheckpointIncrement, float eps,
                                 std::vector<ExperimentResult>& results,
                                 const PointCloud& cloud);
