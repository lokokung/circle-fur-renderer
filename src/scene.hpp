#pragma once

#include <libconfig.h++>

#include <vector>

#include "camera.hpp"
#include "light.hpp"
#include "logger.hpp"
#include "sphere.hpp"

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
