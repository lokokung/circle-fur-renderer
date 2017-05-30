#pragma once

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#define _USE_MATH_DEFINES

#include <Eigen/Eigen>

using namespace Eigen;

struct Vec3i {
    float i, j, k;

    Vec3i(int i, int j, int k);
    Vec3i(int *v);
};

struct Vec3f {
    float x, y, z;

    Vec3f(float x, float y, float z);
    Vec3f(float *v);
};

struct Vec4f {
    float x, y, z, w;

    Vec4f(float x, float y, float z, float w);
    Vec4f(float *v);
};
