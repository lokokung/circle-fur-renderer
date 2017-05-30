#pragma once

#ifdef __CUDA_ARCH__

#include <cuda_runtime.h>
#define CUDA_CALLABLE __host__ __device__

#else

#define CUDA_CALLABLE

#endif
