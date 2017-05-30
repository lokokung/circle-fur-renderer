#pragma once

#include <Eigen/Eigen>

#include "cuda_header.cuh"

class Light {
public:
    // Light position
    Eigen::Vector3f position;
    // Light color
    Eigen::Vector3f color;
    // Attenuation constant
    float k;

    CUDA_CALLABLE Light(float px, float py, float pz,
                        float r, float g, float b, float k);
};
