#pragma once

#include <cstdlib>

#include <Eigen/Eigen>

#include <libnoise/noise.h>

#include "logger.hpp"
#include "sphere.hpp"

#define STEP_SIZE 0.01
#define TILING 512
#define EXPAND 10.0

#define RANDF(min, max) min + ((float) rand() / RAND_MAX) * (max - min)

class FurUtils {
private:
    noise::module::Perlin u_noise, v_noise, l_noise;
    float u_base, v_base, l_base;
    
public:
    FurUtils();
    float fur_density(Sphere* sphere, Eigen::Vector3f point,
                      int tiling, float expand);
    float fur_shadow(Sphere* sphere, Eigen::Vector3f point,
                     Eigen::Vector3f lp,
                     int tiling, float expand);
};
