#pragma once

// Standard lib
#include <cmath>
#include <vector>


namespace Functions {
  // Penalized interpolation function for 0-1 range and its derivative
  // y=x for iPenal= 0
  // Sag halfway below for +3.2
  // Bulge halfway above for -3.2
  // https://www.desmos.com/calculator/ztw49qndpf
  // https://www.wolframalpha.com/input?i=derivative+y%3D+x%2F%281%2Bc-cx%29
  template <std::floating_point element_type>
  inline element_type PenalInterpo(element_type const iVal, element_type const iPenal, bool const iDeriv) {
    const element_type coeff= std::pow(2.0, iPenal) - 1.0;                                             // c= 2^k - 1
    if (iDeriv) return (coeff + 1.0) / ((coeff * iVal - coeff - 1.0) * (coeff * iVal - coeff - 1.0));  // y= (c+1) / (cx-c-1)^2
    else return iVal / (1.0 + coeff - coeff * iVal);                                                   // y= x / (1+c-cx)
  }
}  // namespace Functions
