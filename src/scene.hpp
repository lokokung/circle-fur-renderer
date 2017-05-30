#pragma once

#include <vector>

#include "../cuda/camera.cuh"
#include "../cuda/light.cuh"
#include "../cuda/sphere.cuh"

#ifndef __CUDACC__

#include <libconfig.h++>
#include "logger.hpp"

#endif

class Scene {
public:
    std::vector<Light*> lights;
    Camera* camera;
    Sphere* sphere;

    Scene(Camera* camera, Sphere* sphere);
    Scene(const char* scene_file = "scenes/default.scn");
    ~Scene();
    void add_light(Light* light);

private:
    void construct(Camera* camera, Sphere* sphere);
};
