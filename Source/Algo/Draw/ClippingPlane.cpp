#include "ClippingPlane.hpp"

// Standard lib
#include <array>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"


void ClippingPlane::DrawClipping(const bool iUse, const int iDim, const double iPos, const bool iSide,
                                 const std::array<double, 3> iBoxMin, const std::array<double, 3> iBoxMax) {
  if (iUse) {
    double plane[4];
    if (iDim == 1) {
      plane[0]= -1.0;
      plane[1]= 0.0;
      plane[2]= 0.0;
      plane[3]= iPos * (iBoxMax[0] - iBoxMin[0]) + iBoxMin[0];
    }
    else if (iDim == 2) {
      plane[0]= 0.0;
      plane[1]= -1.0;
      plane[2]= 0.0;
      plane[3]= iPos * (iBoxMax[1] - iBoxMin[1]) + iBoxMin[1];
    }
    else {
      plane[0]= 0.0;
      plane[1]= 0.0;
      plane[2]= -1.0;
      plane[3]= iPos * (iBoxMax[2] - iBoxMin[2]) + iBoxMin[2];
    }

    if (iSide) {
      plane[0]= -plane[0];
      plane[1]= -plane[1];
      plane[2]= -plane[2];
      plane[3]= -plane[3];
    }
    glClipPlane(GL_CLIP_PLANE0, plane);
    glEnable(GL_CLIP_PLANE0);
  }
  else {
    glDisable(GL_CLIP_PLANE0);
  }
}
