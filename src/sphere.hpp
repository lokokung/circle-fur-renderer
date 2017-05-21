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
    // Hair attenuation constant
    float hair_atten;

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
           float phong, float hair_length, float hair_atten);

    // Compute spherical coordinates in radians
    void spherical_coord(Eigen::Vector3f point, float* theta, float* phi);
    
    int intersect_base(Eigen::Vector3f start, Eigen::Vector3f direction,
                       float* t_minus, float* t_plus);
    int intersect_hair(Eigen::Vector3f start, Eigen::Vector3f direction,
                       float* t_minus, float* t_plus);

private:
    int intersect_sphere(Eigen::Vector3f start, Eigen::Vector3f direction,
                         float r, float* t_minus, float* t_plus);
};
