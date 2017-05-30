#pragma once

#include "cuda_header.cuh"

#include "camera.cuh"
#include "light.cuh"
#include "sphere.cuh"
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

private:
    void cuda_raytrace_base(Scene* scene, int xres, int yres,
                            Vec3f *screen, int n_threads);
};
