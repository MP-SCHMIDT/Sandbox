#include "LinearSparseSolverGPU_Kernel.hpp"

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


std::string LinearSparseSolverGPU_Kernel::get_opencl_c_code() {
  std::string r= LinearSparseSolverGPU_Kernel::opencl_c_container();
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
std::string LinearSparseSolverGPU_Kernel::opencl_c_container() {
  return R(

    kernel void kernel_MatrixMult(global const int* arrMatRowStart, global const int* arrMatCol, global const float* arrMatVal,
                                  global const float* arrVec, global float* arrResult, const int N) {
      const unsigned int gid= get_global_id(0);
      if (gid >= N) return;
      float sum= 0.0;
      for (unsigned int idxNNZ= arrMatRowStart[gid]; idxNNZ < arrMatRowStart[gid+1]; idxNNZ++) {
        sum+= arrMatVal[idxNNZ] * arrVec[arrMatCol[idxNNZ]];
      }
      arrResult[gid]= sum;
    }

    kernel void kernel_VecMult(global const float* arrA, global const float* arrB, global float* arrResult, const int N) {
      const unsigned int gid= get_global_id(0);
      if (gid >= N) return;
      arrResult[gid]= arrA[gid] * arrB[gid];
    }

    kernel void kernel_AddScalarMult(global const float* arrA, const float scalar, global const float* arrB, global float* arrResult, const int N) {
      const unsigned int gid= get_global_id(0);
      if (gid >= N) return;
      arrResult[gid]= arrA[gid] + scalar * arrB[gid];
    }

    kernel void kernel_ArrayCopy(global const float* arrA, global float* arrResult, const int N) {
      const unsigned int gid= get_global_id(0);
      if (gid >= N) return;
      arrResult[gid]= arrA[gid];
    }

    kernel void kernel_DotProdPartialReduc(global const float* arrA, global const float* arrB, global float* arrPartial, const int N) {
      local float localSums[64];
      const unsigned int gid= get_global_id(0);
      const unsigned int lid= get_local_id(0);
      const unsigned int group_size= get_local_size(0);
      localSums[lid]= arrA[gid] * arrB[gid];
      for (unsigned int stride= group_size / 2; stride > 0; stride/= 2) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (lid < stride)
          localSums[lid]+= localSums[lid + stride];
      }
      if (lid == 0 && gid < N) // Must check global id here to allow all work items to go through the CLK_LOCAL_MEM_FENCE
        arrPartial[get_group_id(0)]= localSums[0];
    }

  );
}  // End of stringified OpenCL C code
