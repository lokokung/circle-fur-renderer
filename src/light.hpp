#pragma once

#include <Eigen/Eigen>

class Light {
public:
    // Light position
    Eigen::Vector3f position;
    // Light color
    Eigen::Vector3f color;
    // Attenuation constant
    float k;

    Light(float px, float py, float pz,
          float r, float g, float b, float k);
};
