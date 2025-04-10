#pragma once

// OpenCL lib and wrapper
#include "OpenCL_Wrapper/opencl.hpp"


struct
{
  bool isDeviceReady= false;
  bool isMemoryReady= false;
  bool isParamsReady= false;
  int N;
  char lock2D;
  char torusPos;
  char torusFor;
  float timestep;
  float gravity;
  float drag;
  float collision;
  float distDrag;
  float distCollision;
  Device deviceLink;
  Kernel kernelFunc_ComputeForces;
  Kernel kernelFunc_Integrate;
  Memory<cl_float4> arrPos;
  Memory<cl_float4> arrVel;
  Memory<cl_float4> arrFor;
} OCL_NBody;
