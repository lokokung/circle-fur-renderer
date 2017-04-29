#include "raytracer.hpp"

void Raytracer::raytrace_base(Scene *scene, int xres, int yres,
                             const char* output) {
    PNGMaker png(xres, yres);

    // Define maximum color
    Eigen::Vector3f max_color;
    max_color << 1.0, 1.0, 1.0;

    // Get camera information
    float near, aspect, fov, h, w;
    near = scene->camera->near;
    aspect = scene->camera->aspect;
    fov = scene->camera->fov;
    h = 2 * near * tan(fov * M_PI / 360.0);
    w = aspect * h;

    // Make basis vectors
    Eigen::Vector3f e1, e2, e3;
    e1 << 0.0, 0.0, -1.0;
    e2 << 1.0, 0.0, 0.0;
    e3 << 0.0, 1.0, 0.0;

    // Get camera position
    Eigen::Vector3f cp = scene->camera->position;
    // Get position of the sphere
    Eigen::Vector3f sp = scene->sphere->position;
    // Get color vectors of the sphere
    Eigen::Vector3f diffuse, specular, ambient;
    float p = scene->sphere->p;
    diffuse = scene->sphere->diffuse;
    specular = scene->sphere->specular;
    ambient = scene->sphere->ambient;
    // Get position and color of light
    Eigen::Vector3f lp = scene->light->position;
    Eigen::Vector3f lc = scene->light->color;

    // Iterate the pixel grid
    for (int i = 0; i < xres; ++i) {
        for (int j = 0; j < yres; ++j) {
            // Define xi and yj
            float xi = (i - xres / 2.0) * (w / xres);
            float yj = (j - yres / 2.0) * (h / yres);
            // Compute direction of the camera ray
            Eigen::Vector3f direction =
                (near * e1) + (xi * e2) + (yj * e3);
            direction = direction / direction.norm();

            // Compute intersection between sphere and camera ray
            float t_minus, t_plus;
            int res = scene->sphere->intersect_base(cp, direction,
                                                    &t_minus, &t_plus);
            
            // Only continue if we actually intersected with something
            if (res) {
                // Compute point of intersection and normal
                Eigen::Vector3f point = cp + t_minus * direction;
                Eigen::Vector3f normal = point - sp;
                normal = normal / normal.norm();

                // Currently NOT handling shadows or attenuation

                // Get light distance
                Eigen::Vector3f ld = lp - point;
                ld = ld / ld.norm();
                
                // Compute diffuse term for lighting
                Eigen::Vector3f ldiffuse =
                    lc * std::max((float) 0.0, normal.dot(ld));

                // Compute specular term for lighting
                Eigen::Vector3f edld = ld - direction;
                edld = edld / edld.norm();
                Eigen::Vector3f lspecular =
                    lc * std::pow(std::max((float) 0.0, normal.dot(edld)), p);

                // Compute ambient term for lighting
                Eigen::Vector3f lambient = lc;

                // Compute total color
                Eigen::Vector3f color =
                    diffuse.cwiseProduct(ldiffuse) +
                    specular.cwiseProduct(lspecular) +
                    ambient.cwiseProduct(lambient);
                color = color.cwiseMin(max_color);

                // Set the png color
                png.setPixel(i, j, color(0,0), color(1, 0), color(2, 0));   
            }
        }
    }

    // Save the image or report an error
    if (png.saveImage(output)) {
        ERROR("Unable to save the image at {0}", output);
    }
    else {
        INFO("Successfully saved image at {0}", output);
    }
}
