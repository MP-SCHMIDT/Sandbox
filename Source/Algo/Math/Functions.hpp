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

  // Logistic sigmoid function with asymptote flattening [-∞ ; +∞] ⇒ [0 ; 1]
  // https://en.wikipedia.org/wiki/Logistic_function
  template <std::floating_point element_type>
  inline element_type LogisticSigmoid(element_type const iVal, element_type const iBandWidth) {
    if (iBandWidth > 0.0) return 1.0 - 1.0 / (1.0 + std::exp(5.0 * iVal / iBandWidth));  // y= 1-1/(e^(5x/w))
    else return (iVal < 0.0) ? 0.0 : 1.0;
  }

  // Smooth step function [-∞ ; +∞] ⇒ [0 ; 1] with exact match and zero gradient at 0=0 and 1=1
  // https://en.wikipedia.org/wiki/Smoothstep
  template <std::floating_point element_type>
  inline element_type SmoothStep(element_type const iVal) {
    if (iVal > 0.0 && iVal < 1.0) return 3.0 * iVal * iVal - 2.0 * iVal * iVal * iVal;  // y= 3x^2 - 2x^3
    else return (iVal < 0.5) ? 0.0 : 1.0;
  }

  // Smooth min or max of a list of numbers
  // Tends to max with alpha positive and min with alpha negative
  // https://en.wikipedia.org/wiki/Smooth_maximum
  template <std::floating_point element_type>
  inline element_type SmoothMinMax(std::vector<element_type> const iVal, element_type const iAlpha) {
    element_type valMinMax= 0.0, normali= 0.0;
    for (auto val : iVal) {
      valMinMax+= val * std::exp(iAlpha * val);
      normali+= std::exp(iAlpha * val);
    }
    return valMinMax / normali;
  }
}  // namespace Functions
