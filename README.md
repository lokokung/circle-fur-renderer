## Project Details

The goal of the project is to write raytracing software to render images of
spheres covered in hairs. Notably, the choice of rendering spheres was
arbitrary but necessary in order to reduce the possible scope of the project
to be feasible in the given time frame. Note that rendering surfaces with
hairs has been historically a computationally expensive endeavor. As a result,
papers with new rendering algorithms that take away complicated geometry
and replace it with volume spaces to generate the hairs have become important
in the research field. Using the algorithm suggested by Kajiya and Kay in
"Rendering Fur with Three Dimensional Textures," this software approaches the
problem in a similar yet novel fashion.

To distinguish between the algorithm suggested in the paper, namely the
texel generation, we have replaced the fur generation preprocessing step
with a procedural generation of the fur. This is accomplished by using a
modified 2D sinusoidal function mapped onto the base of the sphere to
represent the density of the fur.

Furthermore, using current GPU technologies, the software allows users to
render scenes using CUDA technologies to pass the raytracing work to the
GPU. This has proven to boost performance by between 300-500x compared to
rendering on the CPU.

## Project Library Dependencies

The software provided is only tested in a Linux environment and is dependent
on the following external libraries:

* libconfig++
* libpng
* libnoise
* libtclap
* libeigen

The libraries above can all be installed via package managers on Linux systems.
On Ubunutu systems the following command can be used to do so:

`sudo apt-get install libxxx-dev`

Of course, the source code also requires that CUDA software is installed and
the necessary compiler, nvcc, is available.

## Building Project

The project comes equipped with a Makefile so that compilation should be a
simple make call at the root of the source code. Note that the compiled
binary will be located in the bin folder.

## Project Capabilities and Usage

The software provided can render images of spheres with and without hairs
either on the CPU or the GPU. Options can all be passed in through the
command line interface. As input, the program expects .scn files, and
defaults to the file at scenes/default.scn if not specified. The program
expects either the --CPU or --GPU flag to determine which processor to use.
It also expects either the --base or --hair flag in order to determine what
rendering algorithm to use. For more detailed command line argument information,
use the --help flag provided.

An example call to the program would look as follows:

`bin/furry_circle -i scenes/default.scn -o imgs/ex.png -r 1000 --GPU --hair`

## Current Project Results

Currently, comparing the runtimes between the CPU and GPU, when rendering just
the base sphere at 1000 resolution, the CPU takes roughly 4 seconds while the
GPU only takes 0.2 seconds, a 20x speed up. When rendering with the hairs, the
CPU takes between 2-3 hours while the GPU can complete the image in under 20
seconds, an almost 400x speed up. Notably, the discrepency between the speed up
is caused by the lower computation power needed for just the base rendering
compared to the similar amount of overhead.

The GPU code has yet to be completely optimized and currently scales poorly
(linearly) in runtime with the number of lights. This could easily be fixed,
but due to time constraints was not fixed in this version. 