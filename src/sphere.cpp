#include "sphere.hpp"

Sphere::Sphere(float px, float py, float pz, float radius,
               float diffuse_r, float diffuse_g, float diffuse_b,
               float specular_r, float specular_g, float specular_b,
               float ambient_r, float ambient_g, float ambient_b,
               float phong, float hair_length, float hair_atten) {
    // Set position
    this->position << px, py, pz;
    // Set radius
    this->radius = radius;
    // Set hair length and attenuation
    this->hair_length = hair_length;
    this->hair_atten = hair_atten;

    // Set colors
    this->diffuse << diffuse_r, diffuse_g, diffuse_b;
    this->specular << specular_r, specular_g, specular_b;
    this->ambient << ambient_r, ambient_g, ambient_b;
    this->p = phong;
}

void Sphere::spherical_coord(Eigen::Vector3f point, float* theta, float* phi) {
    // *phi = std::acos(point(2, 0) / radius);
    // *theta = std::atan2(point(1, 0), point(0, 0));

    *phi = std::acos(point(1, 0) / radius);
    *theta = std::atan2(point(2, 0), point(0, 0));
}

int Sphere::intersect_sphere(Eigen::Vector3f start, Eigen::Vector3f direction,
                            float r, float* t_minus, float* t_plus) {
    // Compute individual a, b, and c for quadratic formula
    float a = direction.dot(direction);
    float b = 2.0 * (direction.dot(start) - direction.dot(position));
    float c = position.dot(position) + start.dot(start) -
        2.0 * start.dot(position) - r * r;

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

int Sphere::intersect_base(Eigen::Vector3f start, Eigen::Vector3f direction,
                           float* t_minus, float* t_plus) {
    return intersect_sphere(start, direction, this->radius, t_minus, t_plus);
}

int Sphere::intersect_hair(Eigen::Vector3f start, Eigen::Vector3f direction,
                           float* t_minus, float* t_plus) {
    return intersect_sphere(start, direction,
                            this->radius + this->hair_length,
                            t_minus, t_plus);
}
