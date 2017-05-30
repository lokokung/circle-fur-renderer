#include "noise_perlin.cuh"

Perlin::Perlin() :
    NoiseModule(),
    m_frequency    (DEFAULT_PERLIN_FREQUENCY   ),
    m_lacunarity   (DEFAULT_PERLIN_LACUNARITY  ),
    m_octaveCount  (DEFAULT_PERLIN_OCTAVE_COUNT),
    m_persistence  (DEFAULT_PERLIN_PERSISTENCE ),
    m_seed         (DEFAULT_PERLIN_SEED) {
}

void Perlin::SetFrequency(double frequency) {
    m_frequency = frequency;
}

void Perlin::SetLacunarity(double lacunarity) {
    m_lacunarity = lacunarity;
}

void Perlin::SetOctaveCount(int octaveCount) {
    m_octaveCount = octaveCount;
}

void Perlin::SetPersistence(double persistence) {
    m_persistence = persistence;
}

void Perlin::SetSeed(int seed) {
    m_seed = seed;
}

double Perlin::GetValue(double x, double y, double z) {
    double value = 0.0;
    double signal = 0.0;
    double curPersistence = 1.0;
    double nx, ny, nz;
    int seed;

    x *= m_frequency;
    y *= m_frequency;
    z *= m_frequency;

    for (int curOctave = 0; curOctave < m_octaveCount; curOctave++) {
        
        // Make sure that these floating-point values have the same range as a 32-
        // bit integer so that we can pass them to the coherent-noise functions.
        nx = MakeInt32Range(x);
        ny = MakeInt32Range(y);
        nz = MakeInt32Range(z);

        // Get the coherent-noise value from the input value and add it to the
        // final result.
        seed = (m_seed + curOctave) & 0xffffffff;
        signal = GradientCoherentNoise3D (nx, ny, nz, seed);
        value += signal * curPersistence;

        // Prepare the next octave.
        x *= m_lacunarity;
        y *= m_lacunarity;
        z *= m_lacunarity;
        curPersistence *= m_persistence;
    }

    return value;
}
