#include "logger.hpp"
#include "raytracer.hpp"
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

    Raytracer raytracer;
    raytracer.raytrace_base(scene, 1000, 1000, "test.png");
    
    return 0;
}
