#include "sphere.hpp"

Sphere::Sphere(float px, float py, float pz, float radius,
               float diffuse_r, float diffuse_g, float diffuse_b,
               float specular_r, float specular_g, float specular_b,
               float ambient_r, float ambient_g, float ambient_b,
               float phong, float hair_length) {
    // Set position
    this->position << px, py, pz;
    // Set radius
    this->radius = radius;
    // Set hair length
    this->hair_length = hair_length;

    // Set colors
    this->diffuse << diffuse_r, diffuse_g, diffuse_b;
    this->specular << specular_r, specular_g, specular_b;
    this->ambient << ambient_r, ambient_g, ambient_b;
    this->p = phong;
}

int Sphere::intersect_base(Eigen::Vector3f start, Eigen::Vector3f direction,
                           float* t_minus, float* t_plus) {
    // Compute individual a, b, and c for quadratic formula
    float a = direction.dot(direction);
    float b = 2.0 * (direction.dot(start) - direction.dot(position));
    float c = position.dot(position) + start.dot(start) -
        2.0 * start.dot(position) - radius * radius;

    // Compute discriminant to check if solutions exists
    float discrim = b * b - 4.0 * a * c;
    if (discrim < 0) {
        // If the discriminant is negative, then we have 0 solutions
        return 0;
    }
    else if (discrim == 0) {
        // If the discriminant is 0, then we have exactly 1 solution
        *t_minus = -b / (2.0 * a);
        return 1;
    }
    else {
        // Otherwise there are two solutions
        *t_minus = (-b - sqrt(discrim)) / 2.0 * a;
        *t_plus = (-b + sqrt(discrim)) / 2.0 * a;
        return 2;
    }

    return 0;
}
