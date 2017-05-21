#include "scene.hpp"

Scene::Scene(Camera* camera, Sphere* sphere) {
    construct(camera, sphere);
}

Scene::Scene(const char* scene_file) {
    // Log information
    INFO("Parsing scene file at {0}", scene_file);
    
    // Read configuration file, report error if necessary
    libconfig::Config cfg;
    try {
        cfg.readFile(scene_file);
    }
    catch (const libconfig::FileIOException &fioex) {
        ERROR("I/O error while reading scene file.");
    }
    catch (const libconfig::ParseException &pex) {
        ERROR("Scene file parse error at {0}:{1} - {2}",
              pex.getFile(), pex.getLine(), pex.getError());
    }

    // Prepare scene data to be read
    const libconfig::Setting &root = cfg.getRoot();

    // Read in camera data and make camera
    const libconfig::Setting &camera_data = root["camera"];
    float cam_x, cam_y, cam_z;
    float near, far, fov, aspect;
    camera_data.lookupValue("px", cam_x);
    camera_data.lookupValue("py", cam_y);
    camera_data.lookupValue("pz", cam_z);
    camera_data.lookupValue("near", near);
    camera_data.lookupValue("far", far);
    camera_data.lookupValue("fov", fov);
    camera_data.lookupValue("aspect", aspect);
    Camera *camera = new Camera(cam_x, cam_y, cam_z,
                                near, far, fov, aspect);

    // Read in light data and make light
    const libconfig::Setting &lights_data = root["lights"];
    int n_lights = lights_data.getLength();

    // Read each light
    for (int i = 0; i < n_lights; ++i) {
        // Get the current light
        const libconfig::Setting &light_data = lights_data[i];

        // Read in the parameters for this light
        float light_x, light_y, light_z;
        float light_r, light_g, light_b, light_k;
        light_data.lookupValue("px", light_x);
        light_data.lookupValue("py", light_y);
        light_data.lookupValue("pz", light_z);
        light_data.lookupValue("r", light_r);
        light_data.lookupValue("g", light_g);
        light_data.lookupValue("b", light_b);
        light_data.lookupValue("k", light_k);
        Light *light = new Light(light_x, light_y, light_z,
                                 light_r, light_g, light_b, light_k);

        // Add light to the scene
        lights.push_back(light);
    }
        
    // Read in sphere data and make sphere
    const libconfig::Setting &sphere_data = root["sphere"];
    float sphere_x, sphere_y, sphere_z, radius;
    float hair_length, hair_atten;
    float diffuse_r, diffuse_g, diffuse_b;
    float specular_r, specular_g, specular_b;
    float ambient_r, ambient_g, ambient_b;
    float phong;
    sphere_data.lookupValue("px", sphere_x);
    sphere_data.lookupValue("py", sphere_y);
    sphere_data.lookupValue("pz", sphere_z);
    sphere_data.lookupValue("radius", radius);
    sphere_data.lookupValue("hair_length", hair_length);
    sphere_data.lookupValue("hair_atten", hair_atten);
    sphere_data.lookupValue("rd", diffuse_r);
    sphere_data.lookupValue("gd", diffuse_g);
    sphere_data.lookupValue("bd", diffuse_b);
    sphere_data.lookupValue("rs", specular_r);
    sphere_data.lookupValue("gs", specular_g);
    sphere_data.lookupValue("bs", specular_b);
    sphere_data.lookupValue("ra", ambient_r);
    sphere_data.lookupValue("ga", ambient_g);
    sphere_data.lookupValue("ba", ambient_b);
    sphere_data.lookupValue("p", phong);
    Sphere *sphere = new Sphere(sphere_x, sphere_y, sphere_z, radius,
                                diffuse_r, diffuse_g, diffuse_b,
                                specular_r, specular_g, specular_b,
                                ambient_r, ambient_g, ambient_b,
                                phong, hair_length, hair_atten);
    construct(camera, sphere);
}

Scene::~Scene() {
    delete camera;
    delete sphere;

    for (uint i = 0; i < lights.size(); ++i) {
        delete lights[i];
    }
}

void Scene::construct(Camera* camera, Sphere* sphere) {
    this->camera = camera;
    this->sphere = sphere;
}

void Scene::add_light(Light* light) {
    lights.push_back(light);
}
