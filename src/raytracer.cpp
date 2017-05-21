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

                // Initialize colors for the algorithm
                Eigen::Vector3f ldiffuse, lspecular;
                ldiffuse << 0.0, 0.0, 0.0;
                lspecular << 0.0, 0.0, 0.0;

                // Iterate over every light in the scene
                for (Light* light : scene->lights) {

                    // Get this lights information
                    Eigen::Vector3f lp = light->position;
                    Eigen::Vector3f lc = light->color;
                
                    // Compute point of intersection and normal
                    Eigen::Vector3f point = cp + t_minus * direction;
                    Eigen::Vector3f normal = point - sp;
                    normal = normal / normal.norm();

                    // Currently NOT handling shadows or attenuation

                    // Get light distance
                    Eigen::Vector3f ld = lp - point;
                    ld = ld / ld.norm();
                
                    // Compute diffuse term for lighting
                    ldiffuse +=
                        lc * std::max((float) 0.0, normal.dot(ld));

                    // Compute specular term for lighting
                    Eigen::Vector3f edld = ld - direction;
                    edld = edld / edld.norm();
                    lspecular +=
                        lc * std::pow(std::max((float) 0.0,
                                               normal.dot(edld)), p);
                }
                    
                // Compute total color
                Eigen::Vector3f color =
                    diffuse.cwiseProduct(ldiffuse) +
                    specular.cwiseProduct(lspecular);
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


void Raytracer::raytrace_hair(Scene* scene, int xres, int yres,
                              const char* output) {
    PNGMaker png(xres, yres);
    // Seed random number generator
    srand(time(0));
    // Set up procedural noise generator
    noise::module::Perlin x_noise;
    noise::module::Perlin y_noise;
    x_noise.SetSeed(time(0));
    y_noise.SetSeed(time(0));

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

    // Iterate the pixel grid
    for (int i = 0; i < xres; ++i) {
        INFO("Currently at x: {0}", i);
        for (int j = 0; j < yres; ++j) {
            // Define xi and yj
            float xi = (i - xres / 2.0) * (w / xres);
            float yj = (j - yres / 2.0) * (h / yres);
            // Compute direction of the camera ray
            Eigen::Vector3f direction =
                (near * e1) + (xi * e2) + (yj * e3);
            direction = direction / direction.norm();

            // Compute intersection between sphere and camera ray
            float t_minus_base = std::numeric_limits<float>::max();
            float t_plus_base = std::numeric_limits<float>::max();
            scene->sphere->intersect_base(cp, direction,
                                          &t_minus_base,
                                          &t_plus_base);
            // Compute intersection with hair volume
            float t_minus_hair = std::numeric_limits<float>::max();
            float t_plus_hair = std::numeric_limits<float>::max();
            int res = scene->sphere->intersect_hair(cp, direction,
                                                    &t_minus_hair,
                                                    &t_plus_hair);
            
            // Only continue if we actually intersected with some hair
            if (res > 1) {
                // Compute points of intersection with hair volume
                float t_min = t_minus_hair;
                float t_max = std::min(t_minus_base, t_plus_hair);

                // Initialize colors for the algorithm
                Eigen::Vector3f ldiffuse, lspecular, lambient;
                ldiffuse << 0.0, 0.0, 0.0;
                lspecular << 0.0, 0.0, 0.0;

                // Initialize density summation
                float sum_rho = 0.0;

                // Iterate over every light in the scene
                for (Light* light : scene->lights) {

                    // Get this lights information
                    Eigen::Vector3f lp = light->position;
                    Eigen::Vector3f lc = light->color;
                    float k = light->k;
                
                    // Iterate steps to compute light model for hairs
                    for (float t_low = t_min; t_low < t_max; t_low += STEP_SIZE) {
                        // Compute upper range for segment
                        float t_high = std::min(t_low + (float) STEP_SIZE, t_max);
                    
                        // Get random distance for Monte Carlo integration
                        float t_rand = RANDF(t_low, t_high);

                        // Compute point on intersecting ray
                        Eigen::Vector3f t_point = cp + t_rand * direction;
                        // Compute normal
                        Eigen::Vector3f t_b = t_point - sp;
                        t_b /= t_b.norm();
                        // Compute point on sphere
                        Eigen::Vector3f t_base = t_b * scene->sphere->radius;
                        // Compute vector from base to light and distance
                        Eigen::Vector3f t_l = lp - t_base;
                        float distance = t_l.norm();
                        t_l /= t_l.norm();
                        // Compute vector from base to camera
                        Eigen::Vector3f t_e = cp - t_base;
                        t_e /= t_e.norm();

                        // Apply light distance attenuation
                        Eigen::Vector3f lc_atten =
                            lc / (1.0 + k * distance * distance);
                        lambient = lc_atten;

                        // Compute density at the point
                        float rho = futils.fur_density(scene->sphere, t_base,
                                                       TILING, EXPAND);
                        rho *= (t_high - t_low) / ((float) STEP_SIZE);

                        // Compute second light ray to handle shadowing
                        float shadow = futils.fur_shadow(scene->sphere, t_point,
                                                         lp, TILING, EXPAND);

                        // Compute transparency and increment density
                        float t_t = std::exp(-scene->sphere->hair_atten * sum_rho);
                        sum_rho += rho;
                    
                        // Compute scaled variables with density calculations
                        float t_factors = t_t * rho * shadow;

                        // Compute diffuse term
                        ldiffuse += lc_atten * t_factors * t_b.cross(t_l).norm();

                        // Compute specular term
                        lspecular +=
                            lc_atten * t_factors *
                            std::pow(t_b.dot(t_l) *
                                     t_b.dot(t_e) +
                                     t_b.cross(t_l).norm() *
                                     t_b.cross(t_e).norm(), p);
                    }
                }

                // Compute total color
                Eigen::Vector3f color =
                    diffuse.cwiseProduct(ldiffuse) +
                    specular.cwiseProduct(lspecular);
                    // ambient.cwiseProduct(lambient);
                color = color.cwiseMin(max_color);

                // INFO("RGB: {0}, {1}, {2}", color(0, 0), color(1, 0), color(2, 0));

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


void Raytracer::raytrace_density(Scene *scene, int xres, int yres,
                                 const char* output) {
    PNGMaker png(xres, yres);

    // Set up procedural noise generator
    noise::module::Perlin x_noise;
    noise::module::Perlin y_noise;
    x_noise.SetSeed(time(0));
    y_noise.SetSeed(time(0));

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
    // Eigen::Vector3f sp = scene->sphere->position;
    // Get color vectors of the sphere
    Eigen::Vector3f diffuse, specular, ambient;
    // float p = scene->sphere->p;
    diffuse = scene->sphere->diffuse;
    specular = scene->sphere->specular;
    ambient = scene->sphere->ambient;
    // Get position and color of light
    // Eigen::Vector3f lp = scene->light->position;
    // Eigen::Vector3f lc = scene->light->color;

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

                // Compute spherical coordinates for point
                float theta, phi;
                scene->sphere->spherical_coord(point, &theta, &phi);

                // Generate noise for theta and phi
                // float dtheta = (M_PI / 8.0) * x_noise.GetValue(t_base(0, 0),
                //                                                t_base(1, 0),
                //                                                t_base(2, 0));
                // float dphi = (M_PI / 8.0) * y_noise.GetValue(t_base(0, 0),
                //                                              t_base(1, 0),
                //                                              t_base(2, 0));

                // Compute u and v for mapping
                float u = (TILING * M_PI) * theta / (2.0 * M_PI);
                float v = (TILING * M_PI) * (std::cos(phi) + 1) / 2.0;

                // Compute noise per hair
                float ui = ((int) (u / M_PI)) + 0.01;
                float vi = ((int) (v / M_PI)) + 0.01;
                float du = (M_PI / 4.0) * x_noise.GetValue(ui, vi, 0.1);
                float dv = (M_PI / 4.0) * y_noise.GetValue(ui, vi, 0.1);

                // Compute density at this point
                // float rho =
                //     std::sin(theta * TILING * 2 + dtheta) / 2.0 +
                //     std::sin(phi * TILING + dphi) / 2.0;
                // rho = std::pow(rho, EXPAND);
                float rho =
                    std::sin(u + du) / 2.0 + std::sin(v + dv) / 2.0;
                rho = std::pow(rho, EXPAND);

                // Compute total color
                Eigen::Vector3f color;
                color << rho, rho, rho;

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
