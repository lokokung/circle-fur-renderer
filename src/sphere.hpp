#pragma once

#include <Eigen/Eigen>

class Sphere {
public:
    // Sphere center location
    Eigen::Vector3f position;
    // Sphere radius
    float radius;
    // Length of hair
    float hair_length;

    // Diffuse color
    Eigen::Vector3f diffuse;
    // Specular color
    Eigen::Vector3f specular;
    // Ambient color
    Eigen::Vector3f ambient;

    // Phong exponent
    float p;

    Sphere(float px, float py, float pz, float radius,
           float diffuse_r, float diffuse_g, float diffuse_b,
           float specular_r, float specular_g, float specular_b,
           float ambient_r, float ambient_g, float ambient_b,
           float phong, float hair_length);
    int intersect_base(Eigen::Vector3f start, Eigen::Vector3f direction,
                       float* t_small, float* t_big);
};
