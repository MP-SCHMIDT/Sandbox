#define ENABLE_OPENCL_KERNEL
#ifdef ENABLE_OPENCL_KERNEL

// OpenCL lib and wrapper
#include "OpenCL_Wrapper/stringification.hpp"
#include "OpenCL_Wrapper/utilities.hpp"


// #### Begin of stringified OpenCL C code ####
string opencl_c_container() {
  return R(

      kernel void kernel_AddVec(global const int* A, global const int* B, global int* C) {
        const uint n= get_global_id(0);  // Get index of current work unit
        C[n]= A[n] + B[n];               // Add scalar values together in global memory
      }


      // // https://dournac.org/info/gpu_sum_reduction
      // kernel void kernel_ReducSum(global const int* input, global int* partialSums, int* localSums) {
      //   uint local_id= get_local_id(0);
      //   uint group_size= get_local_size(0);
      //   // Copy from global to local memory
      //   localSums[local_id]= input[get_global_id(0)];
      //   // Loop for computing localSums : divide WorkGroup into 2 parts
      //   for (uint stride= group_size / 2; stride > 0; stride/= 2) {
      //     // Waiting for each 2x2 addition into given workgroup
      //     barrier(CLK_LOCAL_MEM_FENCE);
      //     // Add elements 2 by 2 between local_id and local_id + stride
      //     if (local_id < stride)
      //       localSums[local_id]+= localSums[local_id + stride];
      //   }
      //   // Write result into partialSums[nWorkGroups]
      //   if (local_id == 0)
      //     partialSums[get_group_id(0)]= localSums[0];
      // }

      // https://stackoverflow.com/a/52434838
      kernel void kernel_ReducSum(global int* A, global int* output, local int* target) {
        const size_t globalId= get_global_id(0);
        const size_t localId= get_local_id(0);
        target[localId]= A[globalId];
        barrier(CLK_LOCAL_MEM_FENCE);
        size_t blockSize= get_local_size(0);
        size_t halfBlockSize= blockSize / 2;
        while (halfBlockSize > 0) {
          if (localId < halfBlockSize) {
            target[localId]+= target[localId + halfBlockSize];
            if ((halfBlockSize * 2) < blockSize) {  // uneven block division
              if (localId == 0) {                   // when localID==0
                target[localId]+= target[localId + (blockSize - 1)];
              }
            }
          }
          barrier(CLK_LOCAL_MEM_FENCE);
          blockSize= halfBlockSize;
          halfBlockSize= blockSize / 2;
        }
        if (localId == 0) {
          output[get_group_id(0)]= target[0];
        }
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
}
// #### End of stringified OpenCL C code ####
#endif
