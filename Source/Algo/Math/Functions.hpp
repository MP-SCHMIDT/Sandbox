#pragma once

// Standard lib
#include <vector>


namespace Functions {

  // Penalized interpolation function for 0-1 range and its derivative
  // Equal to y=x for iPenal= 0
  // Sag halfway below for +8
  // Bulge halfway above for -8
  // https://www.desmos.com/calculator/arz1krj34n
  template <std::floating_point element_type>
  inline element_type PenalInterpo(element_type const iVal, element_type const iPenal, bool const iDeriv) {
    if (iPenal >= 0.0) {
      // https://www.wolframalpha.com/input?i=derivative+y%3D+x%2F%281%2Bk-kx%29
      if (iDeriv) return (iPenal + 1.0) / ((iPenal * iVal - iPenal - 1.0) * (iPenal * iVal - iPenal - 1.0));  // (k+1) / (kx-k-1)^2
      else return iVal / (1.0 + iPenal - iPenal * iVal);                                                      // x / (1+k-kx)
    }
    else {
      // https://www.wolframalpha.com/input?i=derivative+y%3D+%28%281-k%29x%29%2F%281-kx%29
      if (iDeriv) return (1.0 - iPenal) / ((iPenal * iVal - 1.0) * (iPenal * iVal - 1.0));  // (1-k) / (kx-1)^2
      else return (1.0 - iPenal) * iVal / (1.0 - iPenal * iVal);                            // (1-k)x / (1-kx)
    }
  }
}  // namespace Functions
