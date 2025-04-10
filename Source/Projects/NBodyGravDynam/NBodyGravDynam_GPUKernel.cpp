#include "NBodyGravDynam_GPUKernel.hpp"

// Standard lib
#include <string>


// Utility function
inline std::string replace(const std::string& s, const std::string& from, const std::string& to) {
	std::string r= s;
	std::size_t p= 0;
	while((p= r.find(from, p)) != std::string::npos) {
		r.replace(p, from.length(), to);
		p+= to.length();
	}
	return r;
}


std::string NBodyGravDynam_GPUKernel::get_opencl_c_code() {
  std::string r= NBodyGravDynam_GPUKernel::opencl_c_container();
  r= replace(r, " ", "\n");                // replace all spaces by new lines
  r= replace(r, "#ifdef\n", "#ifdef ");    // except for the arguments after some preprocessor options that need to be in the same line
  r= replace(r, "#ifndef\n", "#ifndef ");  //
  r= replace(r, "#define\n", "#define ");  // #define with two arguments will not work
	r= replace(r, "#undef\n", "#undef ");    //
  r= replace(r, "#if\n", "#if ");          // don't leave any spaces in arguments
  r= replace(r, "#elif\n", "#elif ");      // don't leave any spaces in arguments
  r= replace(r, "#pragma\n", "#pragma ");
  return "\n" + r;
}


// Evil stringification macro, similar syntax to raw string R"(...)"
// Note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(
#define R(...) std::string(" " #__VA_ARGS__ " ")


// OpenCL lib and wrapper
#include "OpenCL_Wrapper/highlight.hpp"


// Begin of stringified OpenCL C code
std::string NBodyGravDynam_GPUKernel::opencl_c_container() {
  return R(

    kernel void kernel_ComputeForces(const char torusFor, const float gravity, const float drag, const float collision,
                                     const float distDrag, const float distCollision,
                                     global float4* arrPos, global float4* arrVel, global float4* arrFor) {
      const uint N= get_global_size(0); // Get number of bodies
      const uint k0= get_global_id(0);  // Get index of curent body
      // Precomputations
      const float distDragInv= 1.0f / distDrag;
      const float distCollisionInv= 1.0f / distCollision;
      const float distCollisionSqr= distCollision * distCollision;
      // Force integration at predictor step
      const float4 p0= arrPos[k0];
      const float4 v0= arrVel[k0];
      float4 f0= (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      for (uint k1= 0; k1 < N; k1++) {
        // Get body position optionally assuming periodic world
        float4 p1= arrPos[k1];
        if (p0.x == p1.x && p0.y == p1.y && p0.z == p1.z) continue;
        if (torusFor > 0) {
          if      ((p1.x - p0.x) >  0.5f) p1.x= p1.x - 1.0f;
          else if ((p1.x - p0.x) < -0.5f) p1.x= p1.x + 1.0f;
          if      ((p1.y - p0.y) >  0.5f) p1.y= p1.y - 1.0f;
          else if ((p1.y - p0.y) < -0.5f) p1.y= p1.y + 1.0f;
          if      ((p1.z - p0.z) >  0.5f) p1.z= p1.z - 1.0f;
          else if ((p1.z - p0.z) < -0.5f) p1.z= p1.z + 1.0f;
        }
        float4 vec= p1 - p0;
        vec.w= 0.0f;
        const float vecNormSqrInv= 1.0f / dot(vec, vec);
        const float vecNorm= fast_length(vec);
        const float4 vecUnit= vec / vecNorm;

        // Add the gravity force
        f0+= vecUnit * gravity * min(distCollisionSqr * vecNormSqrInv, vecNorm * distCollisionInv);

        // Add the drag force
        if (drag > 0.0f && vecNorm < distDrag)
          f0-= (v0 - arrVel[k1]) * drag * (distDrag - vecNorm) * distDragInv;

        // Add the collision force
        if (collision > 0.0f && vecNorm < distCollision)
          f0-= vecUnit * collision * (distCollision - vecNorm) * distCollisionInv;
      }
      // Save the resulting force
      arrFor[k0]= f0;
    }

    kernel void kernel_Integrate(const char lock2D, const char torusPos, const float timestep,
                                 global float4* arrPos, global float4* arrVel, global float4* arrFor) {
      const uint k0= get_global_id(0);
      float4 v0= arrVel[k0] + timestep * arrFor[k0];
      float4 p0= arrPos[k0] + timestep * v0;
      if (lock2D > 0) {
        p0.x= 0.49f;
        v0.x= 0.0f;
      }
      if (torusPos > 0) {
        if (p0.x < 0.0f || p0.x > 1.0f) p0.x= p0.x - floor(p0.x);
        if (p0.y < 0.0f || p0.y > 1.0f) p0.y= p0.y - floor(p0.y);
        if (p0.z < 0.0f || p0.z > 1.0f) p0.z= p0.z - floor(p0.z);
      }
      arrVel[k0]= v0;
      arrPos[k0]= p0;
    }
  );
}  // End of stringified OpenCL C code
