#include "light.cuh"

Light::Light(float px, float py, float pz,
             float r, float g, float b, float k) {
    // Set light position
    this->position << px, py, pz;

    // Set light color and attenuation
    this->color << r, g, b;
    this->k = k;
}
