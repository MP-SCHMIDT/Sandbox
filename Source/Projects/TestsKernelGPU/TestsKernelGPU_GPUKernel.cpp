#include "TestsKernelGPU_GPUKernel.hpp"

// Standard lib
#include <string>

// OpenCL lib and wrapper
#include "OpenCL_Wrapper/utilities.hpp"

// OpenCL Syntax highlighting
#include "OpenCL_Wrapper/highlight.hpp"


// Evil stringification macro, similar syntax to raw string R"(...)"
// Note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(
#define R(...) std::string(" " #__VA_ARGS__ " ")


std::string TestsKernelGPU_GPUKernel::get_opencl_c_code() {
  std::string r= TestsKernelGPU_GPUKernel::opencl_c_container();
  r= replace(r, " ", "\n");              // replace all spaces by new lines
  r= replace(r, "#ifdef\n", "#ifdef ");  // except for the arguments after some preprocessor options that need to be in the same line
  r= replace(r, "#ifndef\n", "#ifndef ");
  r= replace(r, "#define\n", "#define ");  // #define with two arguments will not work
  r= replace(r, "#if\n", "#if ");          // don't leave any spaces in arguments
  r= replace(r, "#elif\n", "#elif ");      // don't leave any spaces in arguments
  r= replace(r, "#pragma\n", "#pragma ");
  return "\n" + r;
}


// Begin of stringified OpenCL C code
std::string TestsKernelGPU_GPUKernel::opencl_c_container() {
  return R(

      kernel void kernel_AddVec(global const int* A, global const int* B, global int* C) {
        const uint n= get_global_id(0);  // Get index of current work unit
        C[n]= A[n] + B[n];               // Add scalar values together in global memory
      }


      // https://dournac.org/info/gpu_sum_reduction
      kernel void kernel_ReducSum(global const int* input, global int* partialSums) {
        local int localSums[64];  // Creating local with hardcoded size because can't find a way to pass local int* to kernel with opencl wrapper
        uint local_id= get_local_id(0);
        uint group_size= get_local_size(0);
        localSums[local_id]= input[get_global_id(0)];                // Copy from global to local memory
        for (uint stride= group_size / 2; stride > 0; stride/= 2) {  // Loop for computing localSums : divide WorkGroup into 2 parts
          barrier(CLK_LOCAL_MEM_FENCE);                              // Waiting for each 2x2 addition into given workgroup
          if (local_id < stride)                                     // Add elements 2 by 2 between local_id and local_id + stride
            localSums[local_id]+= localSums[local_id + stride];
        }
        if (local_id == 0)  // Write result into partialSums[nWorkGroups]
          partialSums[get_group_id(0)]= localSums[0];
      }


      // https://dournac.org/info/nbody_tutorial
      kernel void kernel_NBodySim(float dt, float eps, float grav,
                                  global float4* pos_old, global float4* pos_new,
                                  global float4* vel_old, global float4* vel_new) {
        const uint N= get_global_size(0);            // Get number of bodies
        const uint k0= get_global_id(0);             // Get index of curent body
        float4 p0= pos_old[k0];                      // Get position of current body
        float4 v= vel_old[k0];                       // Get velocity of current body
        float4 a= (float4)(0.0f, 0.0f, 0.0f, 0.0f);  // Initialize acceleration of current body

        for (uint k1= 0; k1 < N; k1++) {                                                           // Loop through other bodies
          const float4 rVec= pos_old[k1] - p0;                                                     // Get vector from current to other body
          const float rNormInv= rsqrt(rVec.x * rVec.x + rVec.y * rVec.y + rVec.z * rVec.z + eps);  // Compute inverse of Euclidian distance
          a+= grav * rNormInv * rNormInv * rNormInv * rVec;                                        // Add acceleration toward other body
        }
        v+= dt * a;                        // Integrate velocity
        p0+= dt * v + 0.5f * dt * dt * a;  // Integrate position
        vel_new[k0]= v;                    // Save new velocity
        pos_new[k0]= p0;                   // Save new position
      }

  );
}  // End of stringified OpenCL C code
