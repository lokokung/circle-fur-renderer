#include "main.hpp"

int main(int argc, char* argv[]) {
    // Create a console logger to print information
    auto console = spdlog::stdout_color_mt("LOG");

    // Surround argument parsing code in try/catch
    try {
        // Set up argument parser
        TCLAP::CmdLine cmd("Fuzzy sphere procedural generation software",
                           ' ', "1.0");

        // Set up arguments
        TCLAP::ValueArg<std::string> input_file_arg("i", "input_file",
                                                    "Input scene file",
                                                    false,
                                                    "scenes/default.scn",
                                                    "filepath");
        TCLAP::ValueArg<std::string> output_file_arg("o", "output_file",
                                                     "Output png file",
                                                     false,
                                                     "imgs/output.png",
                                                     "filepath");
        
        TCLAP::SwitchArg cpu_arg("", "CPU", "Render using CPU");
        TCLAP::SwitchArg gpu_arg("", "GPU", "Render using GPU");

        TCLAP::SwitchArg base_arg("", "base", "Render image without hairs");
        TCLAP::SwitchArg hair_arg("", "hair", "Render image with hairs");

        TCLAP::ValueArg<int> resolution_arg("r", "resolution",
                                            "Resolution for square image",
                                            false, 1000, "int");

        // Add arguments to parser
        cmd.add(input_file_arg);
        cmd.add(output_file_arg);
        cmd.add(resolution_arg);
        cmd.xorAdd(cpu_arg, gpu_arg);
        cmd.xorAdd(base_arg, hair_arg);

        // Parse the arguments
        cmd.parse(argc, argv);

        // Read arguments
        std::string input = input_file_arg.getValue();
        std::string output = output_file_arg.getValue();
        int resolution = resolution_arg.getValue();

        // Parse scene
        Scene* scene = new Scene(input.c_str());

        // Determine which raytracer to use
        if (cpu_arg.isSet()) {
            Raytracer raytracer;
            // Determine base or hair render
            if (base_arg.isSet()) {
                raytracer.raytrace_base(scene, resolution, resolution,
                                        output.c_str());
            }
            else if (hair_arg.isSet()) {
                raytracer.raytrace_hair(scene, resolution, resolution,
                                        output.c_str());
            }
            else {
                // Should never reach here
                ERROR("Unknown error has occurred");
                exit(-1);
            }
        }
        else if (gpu_arg.isSet()) {
            GPURaytracer raytracer;
            // Determine base or hair render
            if (base_arg.isSet()) {
                raytracer.raytrace_base(scene, resolution, resolution,
                                        output.c_str());
            }
            else if (hair_arg.isSet()) {
                raytracer.raytrace_hair(scene, resolution, resolution,
                                        output.c_str());
            }
            else {
                // Should never reach here
                ERROR("Unknown error has occurred");
                exit(-1);
            }
        }
        else {
            // Should never reach here
            ERROR("Unknown error has occurred");
            exit(-1);
        }
    }
    
    catch (TCLAP::ArgException &e) {
        // If an error is detected log and exit with non-zero error code
        ERROR("TCLAP Error: " + e.error() + " for arg " + e.argId());
        exit(-1);
    }

    // Done, so return 0
    return 0;
}
