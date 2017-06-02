#pragma once

#include <cstdlib>
#include <curand.h>
#include <curand_kernel.h>
#include <float.h>

#include <Eigen/Eigen>

#include "cuda_header.cuh"
#include "noise_perlin.cuh"
#include "sphere.cuh"

#define STEP_SIZE 0.02
#define TILING 256
#define EXPAND 10.0

#define CUDA_RANDF(min, max, st) min + curand_uniform(st) * (max - min)

enum NoiseType {
    U_NOISE,
    V_NOISE,
    N_NOISE_TYPE
};

class GPUFurUtils {
private:
    NoiseModule** modules;
    float *bases;
    curandState *state;
    
public:
    CUDA_CALLABLE GPUFurUtils();
    CUDA_CALLABLE void set_noise_state(curandState* state);
    CUDA_CALLABLE void set_noise_modules(NoiseModule** modules);
    CUDA_CALLABLE void set_noise_bases(float* bases);
    
    __device__ float fur_density(Sphere* sphere, Eigen::Vector3f point,
                                 int tiling, float expand);
    __device__ float fur_shadow(Sphere* sphere, Eigen::Vector3f point,
                                Eigen::Vector3f lp, curandState* state,
                                int tiling, float expand);
};
