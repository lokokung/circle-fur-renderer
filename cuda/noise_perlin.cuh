#pragma once

#include "noise_module_base.cuh"

// Define Perlin Module default parameters
#define DEFAULT_PERLIN_FREQUENCY 1.0
#define DEFAULT_PERLIN_LACUNARITY 2.0
#define DEFAULT_PERLIN_OCTAVE_COUNT 6
#define DEFAULT_PERLIN_PERSISTENCE 0.5
#define DEFAULT_PERLIN_SEED 0
#define PERLIN_MAX_OCTAVE 30

class Perlin: public NoiseModule {
public:
    CUDA_CALLABLE Perlin();
    CUDA_CALLABLE void SetFrequency(double frequency);
    CUDA_CALLABLE void SetLacunarity(double lacunarity);
    CUDA_CALLABLE void SetPersistence(double persistence);
    CUDA_CALLABLE void SetOctaveCount(int octaveCount);
    CUDA_CALLABLE void SetSeed(int seed);

    CUDA_CALLABLE double GetValue(double x, double y, double z);

protected:
    double m_frequency, m_lacunarity, m_persistence;
    int m_octaveCount, m_seed;
};
