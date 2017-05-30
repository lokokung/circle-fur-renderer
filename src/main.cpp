#include "logger.hpp"
#include "raytracer.hpp"
#include "../cuda/gpu_raytracer.cuh"
#include "scene.hpp"

int main(int argc, char* argv[]) {
    // Create a console logger to print information
    auto console = spdlog::stdout_color_mt("LOG");

    // Read scene
    Scene* scene = new Scene();

    // Eigen::Vector3f start;
    // Eigen::Vector3f direction;
    // start << 0.0, 6.0, 10.0;
    // direction << 0.0, 0.0, -1.0;

    // float t1, t2;
    // scene->sphere->intersect_base(start, direction,
    //                               &t1, &t2);
    // INFO("Solutions: {0}, {1}", t1, t2);

    GPURaytracer raytracer;
    raytracer.raytrace_base(scene, 1000, 1000, "imgs/test.png");

    // Eigen::Vector3f point;
    // point << 0.0, -2.0, 0.0;
    // float theta, phi;
    // scene->sphere->spherical_coord(point, &theta, &phi);
    // INFO("Theta: {0}", theta);
    // INFO("Phi: {0}", phi);
    
    return 0;
}
