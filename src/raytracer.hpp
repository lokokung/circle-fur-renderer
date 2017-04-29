#pragma once

#include <algorithm>

#include "logger.hpp"
#include "png_maker.hpp"
#include "scene.hpp"

class Raytracer {
public:
    // Raytracing methods
    void raytrace_base(Scene* scene, int xres, int yres,
                      const char* output);
};
