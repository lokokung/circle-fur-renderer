#include "gpu_raytracer.cuh"

__global__ void raytrace_base_kernel(Camera* camera, Light* light,
                                     Sphere* sphere, int xres, int yres,
                                     Vec3f *screen) {
    // Get current thread id
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    // Verify that we are within the image
    if (index < xres * yres - 1) {

        // Define maximum color
        Eigen::Vector3f max_color;
        max_color << 1.0, 1.0, 1.0;

        // Get camera information
        float near, aspect, fov, h, w;
        near = camera->near;
        aspect = camera->aspect;
        fov = camera->fov;
        h = 2 * near * tanf(fov * M_PI / 360.0);
        w = aspect * h;

        // Make basis vectors
        Eigen::Vector3f e1, e2, e3;
        e1 << 0.0, 0.0, -1.0;
        e2 << 1.0, 0.0, 0.0;
        e3 << 0.0, 1.0, 0.0;

        // Get camera position
        Eigen::Vector3f cp = camera->position;
        // Get position of the sphere
        Eigen::Vector3f sp = sphere->position;
        // Get color vectors of the sphere
        Eigen::Vector3f diffuse, specular;
        float p = sphere->p;
        diffuse = sphere->diffuse;
        specular = sphere->specular;

        // Convert id into x and y in the pixel grid
        int i = index % xres;
        int j = index / xres;

        // Define xi and yj
        float xi = (i - xres / 2.0) * (w / xres);
        float yj = (j - yres / 2.0) * (h / yres);
        // Compute direction of the camera ray
        Eigen::Vector3f direction =
            (near * e1) + (xi * e2) + (yj * e3);
        direction = direction / direction.norm();

        // Compute intersection between sphere and camera ray
        float t_minus, t_plus;
        int res = sphere->intersect_base(cp, direction,
                                         &t_minus, &t_plus);
            
        // Only continue if we actually intersected with something
        if (res) {

            // Initialize colors for the algorithm
            Eigen::Vector3f ldiffuse, lspecular;
            ldiffuse << 0.0, 0.0, 0.0;
            lspecular << 0.0, 0.0, 0.0;

        
            // Get this lights information
            Eigen::Vector3f lp = light->position;
            Eigen::Vector3f lc = light->color;
            float k = light->k;
                
            // Compute point of intersection and normal
            Eigen::Vector3f point = cp + t_minus * direction;
            Eigen::Vector3f normal = point - sp;
            normal = normal / normal.norm();

            // Get light distance
            Eigen::Vector3f ld = lp - point;
            float distance = ld.norm();
            ld = ld / ld.norm();

            // Apply attenuation
            Eigen::Vector3f lc_atten =
                            lc / (1.0 + k * distance * distance);
                
            // Compute diffuse term for lighting
            ldiffuse +=
                lc * fmaxf((float) 0.0, normal.dot(ld));

            // Compute specular term for lighting
            Eigen::Vector3f edld = ld - direction;
            edld = edld / edld.norm();
            lspecular +=
                lc * powf(fmaxf((float) 0.0,
                                normal.dot(edld)), p);
                                   
            // Compute total color
            Eigen::Vector3f color =
                diffuse.cwiseProduct(ldiffuse) +
                specular.cwiseProduct(lspecular);
            color = color.cwiseMin(max_color);
            
            // Set the png color
            screen[index].x += color(0, 0);
            screen[index].y += color(1, 0);
            screen[index].z += color(2, 0);
        }
    }
}
                                     
void GPURaytracer::cuda_raytrace_base(Scene* scene, int xres, int yres,
                                      Vec3f* screen, int n_threads) {
    // Compute kernel parameters
    int n_blocks = ((xres * yres) / n_threads) + 1;
    
    // Declare GPU side data
    Camera* camera;
    Light* light;
    Sphere* sphere;
    Vec3f* dev_screen;
    
    // Allocate GPU memory
    cudaMalloc((void**) &camera, sizeof(Camera));
    cudaMalloc((void**) &light, sizeof(Light));
    cudaMalloc((void**) &sphere, sizeof(Sphere));
    cudaMalloc((void**) &dev_screen, xres * yres * sizeof(Vec3f));

    // Initialize GPU memory
    cudaMemcpy(camera, scene->camera, sizeof(Camera),
               cudaMemcpyDefault);
    cudaMemcpy(sphere, scene->sphere, sizeof(Sphere),
               cudaMemcpyDefault);
    cudaMemset(dev_screen, 0, xres * yres * sizeof(Vec3f));

    // Iterate over each light
    for (Light* l : scene->lights) {
        // Copy light over to GPU
        cudaMemcpy(light, l, sizeof(Light),
                   cudaMemcpyDefault);

        // Call kernel to add this lights effects
        raytrace_base_kernel<<<n_blocks, n_threads>>>(camera, light, sphere,
                                                      xres, yres, dev_screen);
    }

    // Copy screen data back to CPU
    cudaMemcpy(screen, dev_screen, xres * yres * sizeof(Vec3f),
               cudaMemcpyDefault);

    // Clean up GPU memory
    cudaFree(camera);
    cudaFree(light);
    cudaFree(sphere);
    cudaFree(dev_screen);
}


__global__ void raytrace_hair_init(GPUFurUtils** futils,
                                   curandState* states, int seed,
                                   NoiseModule** modules,
                                   int* module_seeds,
                                   float* bases,
                                   int xres, int yres) {

    // Get current thread id
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    // Set up CUDA random library states
    if (index < xres * yres - 1) {
        curand_init(seed, index, 0, &states[index]);
    }
        
    // Use first thread only to set up other state information
    if (index == 0) {
        // Initialize utilities
        *futils = new GPUFurUtils();
        
        // Set bases location and modules location
        (*futils)->set_noise_bases(bases);
        (*futils)->set_noise_modules(modules);
        
        // Initialize noise modules with seeds
        for (int i = 0; i < N_NOISE_TYPE; ++i) {
            modules[i] = new Perlin();
            ((Perlin*) modules[i])->SetSeed(module_seeds[i]);
        }
    }
}

__global__ void raytrace_hair_clean(GPUFurUtils** futils,
                                    NoiseModule** modules) {
    // Free memory used for the utilities
    for (int i = 0; i < N_NOISE_TYPE; ++i) {
        delete modules[i];
    }
    delete (*futils);
}

__global__ void raytrace_hair_kernel(Camera* camera, Light* light,
                                     Sphere* sphere, GPUFurUtils** futils,
                                     curandState* states,
                                     int xres, int yres,
                                     Vec3f *screen) {
    // Get current thread id
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    // Verify that we are within the image
    if (index < xres * yres - 1) {
        
        // Define maximum color
        Eigen::Vector3f max_color;
        max_color << 1.0, 1.0, 1.0;

        // Get camera information
        float near, aspect, fov, h, w;
        near = camera->near;
        aspect = camera->aspect;
        fov = camera->fov;
        h = 2 * near * tanf(fov * M_PI / 360.0);
        w = aspect * h;

        // Make basis vectors
        Eigen::Vector3f e1, e2, e3;
        e1 << 0.0, 0.0, -1.0;
        e2 << 1.0, 0.0, 0.0;
        e3 << 0.0, 1.0, 0.0;

        // Get camera position
        Eigen::Vector3f cp = camera->position;
        // Get position of the sphere
        Eigen::Vector3f sp = sphere->position;
        // Get color vectors of the sphere
        Eigen::Vector3f diffuse, specular;
        float p = sphere->p;
        diffuse = sphere->diffuse;
        specular = sphere->specular;

        // Convert id into x and y in the pixel grid
        int i = index % xres;
        int j = index / xres;

        // Define xi and yj
        float xi = (i - xres / 2.0) * (w / xres);
        float yj = (j - yres / 2.0) * (h / yres);
        // Compute direction of the camera ray
        Eigen::Vector3f direction =
            (near * e1) + (xi * e2) + (yj * e3);
        direction = direction / direction.norm();

        // Compute intersection between sphere and camera ray
        float t_minus_base = FLT_MAX;
        float t_plus_base = FLT_MAX;
        sphere->intersect_base(cp, direction,
                               &t_minus_base,
                               &t_plus_base);
        // Compute intersection with hair volume
        float t_minus_hair = FLT_MAX;
        float t_plus_hair = FLT_MAX;
        int res = sphere->intersect_hair(cp, direction,
                                         &t_minus_hair,
                                         &t_plus_hair);
            
        // Only continue if we actually intersected with some hair
        if (res > 1) {
            // Compute points of intersection with hair volume
            float t_min = t_minus_hair;
            float t_max = min(t_minus_base, t_plus_hair);
            
            // Initialize colors for the algorithm
            Eigen::Vector3f ldiffuse, lspecular;
            ldiffuse << 0.0, 0.0, 0.0;
            lspecular << 0.0, 0.0, 0.0;
            
            // Initialize density summation
            float sum_rho = 0.0;

            // Get this lights information
            Eigen::Vector3f lp = light->position;
            Eigen::Vector3f lc = light->color;
            float k = light->k;
                
            // Iterate steps to compute light model for hairs
            for (float t_low = t_min; t_low < t_max; t_low += STEP_SIZE) {
                // Compute upper range for segment
                float t_high = min(t_low + (float) STEP_SIZE, t_max);
                
                // Get random distance for Monte Carlo integration
                float t_rand = CUDA_RANDF(t_low, t_high, states + index);

                // Compute point on intersecting ray
                Eigen::Vector3f t_point = cp + t_rand * direction;
                // Compute normal
                Eigen::Vector3f t_b = t_point - sp;
                t_b /= t_b.norm();
                // Compute point on sphere
                Eigen::Vector3f t_base = t_b * sphere->radius;
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

                // Compute density at the point
                float rho = (*futils)->fur_density(sphere, t_base,
                                                   TILING, EXPAND);
                rho *= (t_high - t_low) / ((float) STEP_SIZE);

                // Compute second light ray to handle shadowing
                float shadow = (*futils)->fur_shadow(sphere, t_point,
                                                     lp, &states[index],
                                                     TILING, EXPAND);

                // Compute transparency and increment density
                float t_t = expf(-sphere->hair_atten * sum_rho);
                sum_rho += rho;
                    
                // Compute scaled variables with density calculations
                float t_factors = t_t * rho * shadow;

                // Compute diffuse term
                ldiffuse += lc_atten * t_factors * t_b.cross(t_l).norm();

                // Compute specular term
                lspecular +=
                    lc_atten * t_factors *
                    powf(t_b.dot(t_l) *
                         t_b.dot(t_e) +
                         t_b.cross(t_l).norm() *
                         t_b.cross(t_e).norm(), p);
            }

            // Compute total color
            Eigen::Vector3f color =
                diffuse.cwiseProduct(ldiffuse) +
                specular.cwiseProduct(lspecular);
            color = color.cwiseMin(max_color);

            // Set the png color
            screen[index].x += color(0, 0);
            screen[index].y += color(1, 0);
            screen[index].z += color(2, 0);
        }
    }
}

void GPURaytracer::cuda_raytrace_hair(Scene* scene, int xres, int yres,
                                      Vec3f* screen, int n_threads) {
    // Compute kernel parameters
    int n_blocks = ((xres * yres) / n_threads) + 1;

    // Declare GPU side data
    Camera* camera;
    Light* light;
    Sphere* sphere;
    Vec3f* dev_screen;
    GPUFurUtils** futils;
    curandState* states;
    NoiseModule** modules;
    int* module_seeds;
    float* bases;
    
    // Allocate GPU memory
    cudaMalloc((void**) &camera, sizeof(Camera));
    cudaMalloc((void**) &light, sizeof(Light));
    cudaMalloc((void**) &sphere, sizeof(Sphere));
    cudaMalloc((void**) &dev_screen, xres * yres * sizeof(Vec3f));
    cudaMalloc((void**) &futils, sizeof(GPUFurUtils*));
    cudaMalloc((void**) &states, xres * yres * sizeof(curandState));
    cudaMalloc((void**) &modules, N_NOISE_TYPE * sizeof(NoiseModule*));
    cudaMalloc((void**) &module_seeds, N_NOISE_TYPE * sizeof(int));
    cudaMalloc((void**) &bases, N_NOISE_TYPE * sizeof(float));

    // Initialize module seeds, and bases (temporarily hard-coded)
    int host_module_seeds[2] = {1, 2};
    cudaMemcpy(module_seeds, host_module_seeds, N_NOISE_TYPE * sizeof(int),
               cudaMemcpyDefault);
    float host_bases[2] = {0.123, 0.672};
    cudaMemcpy(bases, host_bases, N_NOISE_TYPE * sizeof(float),
               cudaMemcpyDefault);

    // Initialize scene information
    cudaMemcpy(camera, scene->camera, sizeof(Camera),
               cudaMemcpyDefault);
    cudaMemcpy(sphere, scene->sphere, sizeof(Sphere),
               cudaMemcpyDefault);
    cudaMemset(dev_screen, 0, xres * yres * sizeof(Vec3f));

    // Initialize utilities and states (current seed hard-coded)
    raytrace_hair_init<<<n_blocks, n_threads>>>(futils, states, 1234,
                                                modules, module_seeds,
                                                bases, xres, yres);

    // Iterate over each light
    for (Light* l : scene->lights) {
        // Copy light over to GPU
        cudaMemcpy(light, l, sizeof(Light),
                   cudaMemcpyDefault);

        // Call kernel to add this lights effects
        raytrace_hair_kernel<<<n_blocks, n_threads>>>(camera, light, sphere,
                                                      futils, states,
                                                      xres, yres, dev_screen);
    }

    // Copy screen data back to CPU
    cudaMemcpy(screen, dev_screen, xres * yres * sizeof(Vec3f),
               cudaMemcpyDefault);

    // Clean up GPU memory
    raytrace_hair_clean<<<1, 1>>>(futils, modules);
    cudaFree(camera);
    cudaFree(light);
    cudaFree(sphere);
    cudaFree(dev_screen);
    cudaFree(futils);
    cudaFree(states);
    cudaFree(modules);
    cudaFree(module_seeds);
    cudaFree(bases);
}
