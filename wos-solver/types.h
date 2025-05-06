#pragma once

// To compile: g++ -std=c++17 -O3 -pedantic -Wall -I./include wos_fileread.cpp -o wos
#include <array>
#include <complex>
#include <vector>
#include "nanoflann.hpp"
// #include <algorithm>
// #include <functional>
// #include <iostream>
// #include <random>
// #include <fstream>


// --- Hashing and equality for complex<float> ---

namespace std {
    template<>
    struct hash<std::complex<double>> {
        std::size_t operator()(const std::complex<double>& v) const {
            auto h1 = std::hash<double>{}(v.real());
            auto h2 = std::hash<double>{}(v.imag());
            return h1 ^ (h2 << 1);  // Combine the hashes
        }
    };
 
    template<>
    struct equal_to<std::complex<double>> {
        bool operator()(const std::complex<double>& lhs, const std::complex<double>& rhs) const {
            return std::abs(lhs.real() - rhs.real()) < 1e-6 &&
                   std::abs(lhs.imag() - rhs.imag()) < 1e-6;
        }
    };
}

// --- Common type aliases ---

// use std::complex to implement 2D vectors
// using Vec2D = std::complex<float>;
using Vec2D = std::complex<double>;

// a segment is just a pair of points
using Segment = std::array<Vec2D, 2>;

// boundary geometry is represented by polylines
using Polyline = std::vector<Vec2D>;

// --- Point cloud structure for KDTree ---

struct PointCloud {
    struct Point {
       double x, y, temperature;
    };
 
    std::vector<Point> pts;
 
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

// Declare a global pointer for the KDTree
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 2> KDTree;


// --- ExperimentResult struct ---

struct ExperimentResult {
    double epsilon;
    int nWalks;
    double l1_error;
    double cumulative_time;
};