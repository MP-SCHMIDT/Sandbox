#include "TestsKernelGPU_kernel.hpp"
// Note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(

// #### Begin of OpenCL C code ####
string opencl_c_container() {
  return R(

      kernel void kernel_Add(global float* A, global float* B, global float* C) {
        const uint n= get_global_id(0);  // Get index of current work unit
        C[n]= A[n] + B[n];               // Add scalar values togeter
      }


      kernel void kernel_NBody(int N, float dt, float eps, float grav,
                               global float4* pos_old, global float4* pos_new,
                               global float4* vel_old, global float4* vel_new) {
        const int k0= get_global_id(0);              // Get index of curent body
        float4 p0= pos_old[k0];                      // Get position of current body
        float4 v= vel_old[k0];                       // Get velocity of current body
        float4 a= (float4)(0.0f, 0.0f, 0.0f, 0.0f);  // Initialize acceleration of current body

        for (int k1= 0; k1 < N; k1++) {                                                      // Loop through other bodies
          float4 p1= pos_old[k1];                                                            // Get position of other body
          float4 rVec= p1 - p0;                                                              // Get vector from current to other body
          float rNormInv= rsqrt(rVec.x * rVec.x + rVec.y * rVec.y + rVec.z * rVec.z + eps);  // Compute inverse of distance
          a+= grav * rNormInv * rNormInv * rNormInv * rVec;                                  // Add acceleration toward other body
        }
        v+= dt * a;                        // Integrate velocity
        p0+= dt * v + 0.5f * dt * dt * a;  // Integrate position
        vel_new[k0]= v;                    // Save new velocity
        pos_new[k0]= p0;                   // Save new position
      }

  );
}
// #### End of OpenCL C code ####
