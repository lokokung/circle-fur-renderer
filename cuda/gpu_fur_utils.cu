#include "gpu_fur_utils.cuh"

GPUFurUtils::GPUFurUtils() {
    state = NULL;
}

void GPUFurUtils::set_noise_state(curandState* state) {
    this->state = state;
}

void GPUFurUtils::set_noise_modules(NoiseModule** modules) {
    this->modules = modules;
}

void GPUFurUtils::set_noise_bases(float* bases) {
    this->bases = bases;
}

__device__ float GPUFurUtils::fur_density(Sphere* sphere, Eigen::Vector3f point,
                                           int tiling, float expand) {
    // Compute spherical coordinates
    float theta, phi;
    sphere->spherical_coord(point, &theta, &phi);

    // Compute u and v for mapping
    float u = (tiling * M_PI) * theta / (2.0 * M_PI);
    float v = (tiling * M_PI) * (cosf(phi) + 1) / 2.0;

    // Compute noise per hair
    float ui = ((int) (u / M_PI)) + bases[U_NOISE];
    float vi = ((int) (v / M_PI)) + bases[V_NOISE];
    float du = (M_PI / 4.0) *
        ((Perlin*) modules[U_NOISE])->GetValue(ui, vi, bases[U_NOISE]);
    float dv = (M_PI / 4.0) *
        ((Perlin*) modules[V_NOISE])->GetValue(ui, vi, bases[V_NOISE]);
    
    // Compute density at this point
    float rho =
        sinf(u + du) / 2.0 + sinf(v + dv) / 2.0;
    rho = powf(rho, expand);

    return rho;
}

__device__ float GPUFurUtils::fur_shadow(Sphere* sphere, Eigen::Vector3f point,
                                         Eigen::Vector3f lp, curandState* state,
                                         int tiling, float expand) {
    // Compute t for the actual point and direction
    Eigen::Vector3f direction = point - lp;
    float t = direction.norm();
    direction /= t;

    // Intersect ray with sphere
    float t_minus_base = FLT_MAX;
    float t_plus_base  = FLT_MAX;
    sphere->intersect_base(lp, direction,
                           &t_minus_base,
                           &t_plus_base);
    // Intersect ray with hairs
    float t_minus_hair = FLT_MAX;
    float t_plus_hair = FLT_MAX;
    int res = sphere->intersect_hair(lp, direction,
                                     &t_minus_hair,
                                     &t_plus_hair);

    // Check that we intersected (should be trivially true) and then
    // iterate across the ray to get sum of densities
    float sum_rho = 0.0;
    if (res > 1 && t_minus_hair < t && t_minus_base > t) {
        // Iterate across the ray
        for (float t_low = t_minus_hair; t_low < t; t_low += STEP_SIZE) {
            // Compute upper range for segment
            float t_high = min(t_low + (float) STEP_SIZE, t);

            // Get random distance for Monte Carlo integration
            float t_rand = CUDA_RANDF(t_low, t_high, state);

            // Compute point on intersecting ray
            Eigen::Vector3f t_point = lp + t_rand * direction;
            // Compute normal
            Eigen::Vector3f t_b = t_point - sphere->position;
            t_b /= t_b.norm();
            // Compute point on sphere
            Eigen::Vector3f t_base = t_b * sphere->radius;

            // Compute density at the point
            float rho = fur_density(sphere, t_base,
                                    TILING, EXPAND);
            rho *= (t_high - t_low) / ((float) STEP_SIZE);

            // Add to total density
            sum_rho += rho;
        }
        
        // Compute transparency factor
        float transparency = expf(-sphere->hair_atten * sum_rho);
    
        return transparency;
    }

    // Otherwise, the light is completely blocked so return a factor of 0
    return 0.0;
}
