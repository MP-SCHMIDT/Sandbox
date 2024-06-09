#pragma once

// #include "TestsKernelGPU.hpp"

// Standard lib
#include <vector>

// OpenCL lib
#include "OpenCL_Wrapper/opencl.hpp"

// Global headers
// #include "Data.hpp"


Device device;
Memory<cl_float4> PosOld;
Memory<cl_float4> PosNew;
Memory<cl_float4> VelOld;
Memory<cl_float4> VelNew;
Kernel kernel_NBody;
