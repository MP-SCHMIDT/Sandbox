#include "TestsKernelGPU_kernel.hpp"
// Note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(

// #### Begin of OpenCL C code ####
string opencl_c_container() { return R(

// Equivalent to "for(uint n=0u; n<N; n++) {", but executed in parallel
kernel void add_kernel(global float* A, global float* B, global float* C) {
	const uint n = get_global_id(0);
	C[n] = A[n]+B[n];
}

);}
// #### End of OpenCL C code ####