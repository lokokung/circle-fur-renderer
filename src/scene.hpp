#pragma once

#include <libconfig.h++>

#include "camera.hpp"
#include "light.hpp"
#include "logger.hpp"
#include "sphere.hpp"

class Scene {
public:
    Light* light;
    Camera* camera;
    Sphere* sphere;

    Scene(Light* light, Camera* camera, Sphere* sphere);
    Scene(const char* scene_file = "scenes/default.scn");
    ~Scene();

private:
    void construct(Light* light, Camera* cameran, Sphere* sphere);
};
