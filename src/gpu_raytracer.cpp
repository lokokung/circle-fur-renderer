#include "../cuda/gpu_raytracer.cuh"

void GPURaytracer::raytrace_base(Scene* scene, int xres, int yres,
                                 const char* output) {
    // Initialize PNG maker
    PNGMaker png(xres, yres);

    // Call GPU helper to pass work to GPU
    cuda_raytrace_base(scene, xres, yres, png.screen, 128);

    // Save the image or report an error
    if (png.saveImage(output)) {
        ERROR("Unable to save the image at {0}", output);
    }
    else {
        INFO("Successfully saved image at {0}", output);
    }
}

void GPURaytracer::raytrace_hair(Scene* scene, int xres, int yres,
                                 const char* output) {
    // Initialize PNG maker
    PNGMaker png(xres, yres);

    // Call GPU helper to pass work to GPU
    cuda_raytrace_hair(scene, xres, yres, png.screen, 128);

    // Save the image or report an error
    if (png.saveImage(output)) {
        ERROR("Unable to save the image at {0}", output);
    }
    else {
        INFO("Successfully saved image at {0}", output);
    }
}
