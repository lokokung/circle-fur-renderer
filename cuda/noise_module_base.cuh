#pragma once

#include "cuda_header.cuh"

/* Noise module class based on libnoise library code. Adapted for
 * CUDA usage. */
class NoiseModule {
public:
    CUDA_CALLABLE NoiseModule();
    CUDA_CALLABLE double GetValue(double x, double y, double z);

protected:
    // CUDA implementation of constant table lookup
    CUDA_CALLABLE double g_randomVectors(int i);
    
    // Functions originally from noisegen from libnoise
    CUDA_CALLABLE double MakeInt32Range(double n);
    CUDA_CALLABLE double GradientNoise3D(double fx, double fy, double fz,
                                         int ix, int iy, int iz,
                                         int seed);
    CUDA_CALLABLE double GradientCoherentNoise3D(double x, double y, double z,
                                                 int seed);

    // Functions originally from interp from libnoise
    CUDA_CALLABLE double LinearInterp(double n0, double n1, double a);
    CUDA_CALLABLE double SCurve3(double a);
    
};
