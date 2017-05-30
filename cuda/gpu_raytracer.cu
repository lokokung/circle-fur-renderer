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
        h = 2 * near * tan(fov * M_PI / 360.0);
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
