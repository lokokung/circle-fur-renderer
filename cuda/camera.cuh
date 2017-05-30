#pragma once

#include <Eigen/Eigen>

#include "cuda_header.cuh"

class Camera {
public:
    // Camera position
    Eigen::Vector3f position;
    // Frustrum information
    float near, far, fov, aspect;

    CUDA_CALLABLE Camera(float px, float py, float pz,
                         float near, float far,
                         float fov, float aspect);
};
