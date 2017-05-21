#include "fur_utils.hpp"

FurUtils::FurUtils() {
    // Seed procedural generators
    srand(time(0));
    u_noise.SetSeed(time(0));
    v_noise.SetSeed(time(0));
    l_noise.SetSeed(time(0));

    // Initialize bases
    u_base = ((float) rand() / RAND_MAX);
    v_base = ((float) rand() / RAND_MAX);
    l_base = ((float) rand() / RAND_MAX);
}

float FurUtils::fur_density(Sphere* sphere, Eigen::Vector3f point,
                            int tiling, float expand) {
    // Compute spherical coordinates
    float theta, phi;
    sphere->spherical_coord(point, &theta, &phi);

    // Compute u and v for mapping
    float u = (tiling * M_PI) * theta / (2.0 * M_PI);
    float v = (tiling * M_PI) * (std::cos(phi) + 1) / 2.0;

    // Compute noise per hair
    float ui = ((int) (u / M_PI)) + u_base;
    float vi = ((int) (v / M_PI)) + v_base;
    float du = (M_PI / 4.0) * u_noise.GetValue(ui, vi, u_base);
    float dv = (M_PI / 4.0) * v_noise.GetValue(ui, vi, v_base);

    // Compute density at this point
    float rho =
        std::sin(u + du) / 2.0 + std::sin(v + dv) / 2.0;
    rho = std::pow(rho, expand);

    return rho;
}

float FurUtils::fur_shadow(Sphere* sphere, Eigen::Vector3f point,
                           Eigen::Vector3f lp,
                           int tiling, float expand) {
    // Compute t for the actual point and direction
    Eigen::Vector3f direction = point - lp;
    float t = direction.norm();
    direction /= t;

    // Intersect ray with sphere
    float t_minus_base = std::numeric_limits<float>::max();
    float t_plus_base = std::numeric_limits<float>::max();
    sphere->intersect_base(lp, direction,
                           &t_minus_base,
                           &t_plus_base);
    // Intersect ray with hairs
    float t_minus_hair = std::numeric_limits<float>::max();
    float t_plus_hair = std::numeric_limits<float>::max();
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
            float t_high = std::min(t_low + (float) STEP_SIZE, t);

            // Get random distance for Monte Carlo integration
            float t_rand = RANDF(t_low, t_high);

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
        float transparency = std::exp(-sphere->hair_atten * sum_rho);
    
        return transparency;
    }

    // Otherwise, the light is completely blocked so return a factor of 0
    return 0.0;
}
