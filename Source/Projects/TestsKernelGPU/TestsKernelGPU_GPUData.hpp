#pragma once

// OpenCL lib and wrapper
#include "OpenCL_Wrapper/opencl.hpp"


struct
{
  int N;
  Device deviceLink;
  Kernel kernelFuncABC;
  Kernel kernelFuncBCA;
  Memory<int> A;
  Memory<int> B;
  Memory<int> C;
} OCL_AddVec;

struct
{
  int N;
  int size;
  Device deviceLink;
  Kernel kernelFunc;
  Memory<int> arr;
  Memory<int> par;
  Memory<int> loc;
} OCL_ReducSum;

struct
{
  bool isDeviceReady= false;
  bool isMemoryReady= false;
  bool isParamsReady= false;
  int N;
  float timestep;
  float epsilon;
  float gravity;
  Device deviceLink;
  Kernel kernelFunc;
  Memory<cl_float4> PosOld;
  Memory<cl_float4> PosNew;
  Memory<cl_float4> VelOld;
  Memory<cl_float4> VelNew;
} OCL_NBody;
