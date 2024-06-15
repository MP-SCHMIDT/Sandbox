#pragma once

// OpenCL lib and wrapper
#include "OpenCL_Wrapper/utilities.hpp"

#define R(...) string(" "#__VA_ARGS__" ") // evil stringification macro, similar syntax to raw string R"(...)"
// Note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(

string opencl_c_container(); // outsourced to kernel.cpp
string get_opencl_c_code() {
	string r = opencl_c_container();
	r = replace(r, " ", "\n"); // replace all spaces by new lines
	r = replace(r, "#ifdef\n", "#ifdef "); // except for the arguments after some preprocessor options that need to be in the same line
	r = replace(r, "#ifndef\n", "#ifndef ");
	r = replace(r, "#define\n", "#define "); // #define with two arguments will not work
	r = replace(r, "#if\n", "#if "); // don't leave any spaces in arguments
	r = replace(r, "#elif\n", "#elif "); // don't leave any spaces in arguments
	r = replace(r, "#pragma\n", "#pragma ");
	return "\n"+r;
}

// OpenCL lib and wrapper
#include "OpenCL_Wrapper/highlight.hpp" // Include after defining the stringification macro
