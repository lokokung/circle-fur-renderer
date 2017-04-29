#include "camera.hpp"

Camera::Camera(float px, float py, float pz,
               float near, float far,
               float fov, float aspect) {
    this->position << px, py, pz;
    
    this->near = near;
    this->far = far;
    this->fov = fov;
    this->aspect = aspect;
}

