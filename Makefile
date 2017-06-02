# Linux Makefile

# Product Names
EXE_NAME = bin/furry_circle
CUDA_OBJ = obj/device.o

# Necessary directories
DIR_PATHS = bin obj

# ------------------------------------------------------------------------------

# CUDA Compiler and Flags
NVCC = nvcc
CUDA_FLAGS = -g -m64 -dc -Wno-deprecated-gpu-targets --std=c++11
CUDA_INCLUDE = -I /usr/local/cuda/include -I /usr/include -I /usr/include/eigen3
CUDA_LIBS = 
CUDA_GENCODES = -gencode arch=compute_20,code=sm_20 \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-gencode arch=compute_50,code=sm_50 \
		-gencode arch=compute_52,code=sm_52 \
		-gencode arch=compute_60,code=sm_60 \
		-gencode arch=compute_60,code=compute_60

# CUDA Source Files
CUDA_FILES = $(wildcard cuda/*.cu)

# CUDA Object Files
CUDA_OBJ_FILES = $(addprefix obj/, $(notdir $(addsuffix .o, $(CUDA_FILES))))

# ------------------------------------------------------------------------------

# CUDA Linker and Flags
CUDA_LINK_FLAGS = -dlink -Wno-deprecated-gpu-targets

# ------------------------------------------------------------------------------

# C++ Compiler and Flags
GPP = g++
FLAGS = -g -Wall -D_REENTRANT -std=c++0x -pthread
INCLUDE = -isystem include -I /usr/include -I /usr/include/eigen3 \
	  -I /usr/local/cuda/include
LIBS = -lconfig++ -lpng -lnoise -L/usr/local/cuda/lib64 -lcudart

# C++ Source Files
CPP_FILES = $(wildcard src/*.cpp)

# C++ Object Files
OBJ_FILES = $(addprefix obj/, $(notdir $(addsuffix .o, $(CPP_FILES))))


# ------------------------------------------------------------------------------
# Make Rules
# ------------------------------------------------------------------------------

all: $(EXE_NAME)

# Make executable
$(EXE_NAME): $(OBJ_FILES) $(CUDA_OBJ) $(CUDA_OBJ_FILES)
	$(GPP) $(FLAGS) -o $@ $(INCLUDE) $^ $(LIBS) 

# Make linked device code
$(CUDA_OBJ): $(CUDA_OBJ_FILES)
	$(NVCC) $(CUDA_LINK_FLAGS) $(CUDA_GENCODES) -o $@ $(CUDA_INCLUDE) $^

# Make output directories
dir:
	mkdir -p $(DIR_PATHS)

# Compile C++ Source Files
obj/%.cpp.o: src/%.cpp | dir
	$(GPP) $(FLAGS) -c -o $@ $(INCLUDE) $< 

# Compile CUDA Source Files
obj/%.cu.o: cuda/%.cu | dir
	$(NVCC) $(CUDA_FLAGS) $(CUDA_GENCODES) -c -o $@ $(CUDA_INCLUDE) $<

src: $(OBJ_FILES)

cuda: $(CUDA_OBJ_FILES) $(CUDA_OBJ)

# Clean everything including temporary Emacs files
clean:
	rm -rf bin/* obj/* $(EXE_NAME)
	rm -f src/*~ scenes/*~ cuda/*~
	rm -f *~

.PHONY: clean
