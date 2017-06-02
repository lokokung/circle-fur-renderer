#pragma once

#include <iostream>

#include "cuda_header.cuh"

#include "camera.cuh"
#include "light.cuh"
#include "sphere.cuh"
#include "gpu_fur_utils.cuh"

#include "noise_perlin.cuh"

#include "../src/scene.hpp"
#include "../src/utilities.hpp"

#ifndef __CUDACC__
#include "../src/logger.hpp"
#include "../src/png_maker.hpp"
#endif

class GPURaytracer {
public:
    // Raytracing methods
    void raytrace_base(Scene* scene, int xres, int yres,
                       const char* output);
    void raytrace_hair(Scene* scene, int xres, int yres,
                       const char* output);

private:
    void cuda_raytrace_base(Scene* scene, int xres, int yres,
                            Vec3f* screen, int n_threads);
    void cuda_raytrace_hair(Scene* scene, int xres, int yres,
                            Vec3f* screen, int n_threads);
};
