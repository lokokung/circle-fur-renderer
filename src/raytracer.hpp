#pragma once

#include <libnoise/noise.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "fur_utils.hpp"
#include "logger.hpp"
#include "png_maker.hpp"
#include "scene.hpp"

class Raytracer {
private:
    FurUtils futils;
public:
    // Raytracing methods
    void raytrace_base(Scene* scene, int xres, int yres,
                      const char* output);
    void raytrace_hair(Scene* scene, int xres, int yres,
                       const char* output);
    void raytrace_density(Scene* scene, int xres, int yres,
                          const char* output);
};
