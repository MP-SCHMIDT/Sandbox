#include "DrawShape.hpp"

// GLUT lib
#include "GL/freeglut.h"


void DrawShape::DrawBox(const float begX, const float begY, const float begZ,
                        const float sizX, const float sizY, const float sizZ,
                        bool const isSolid) {
  glPushMatrix();
  glTranslatef(begX, begY, begZ);
  glScalef(sizX, sizY, sizZ);
  glTranslatef(0.5f, 0.5f, 0.5f);
  if (isSolid) glutSolidCube(1.0);
  else glutWireCube(1.0);
  glPopMatrix();
}
