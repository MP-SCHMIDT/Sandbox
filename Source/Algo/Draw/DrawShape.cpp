#include "DrawShape.hpp"

// GLUT lib
#include "freeglut/include/GL/freeglut.h"


void DrawShape::DrawBoxPosPos(const float begX, const float begY, const float begZ,
                              const float endX, const float endY, const float endZ, bool const isSolid) {
  glPushMatrix();
  glTranslatef(begX, begY, begZ);
  glScalef(endX - begX, endY - begY, endZ - begZ);
  glTranslatef(0.5f, 0.5f, 0.5f);
  if (isSolid) glutSolidCube(1.0);
  else glutWireCube(1.0);
  glPopMatrix();
}


void DrawShape::DrawBoxPosSiz(const float begX, const float begY, const float begZ,
                              const float sizX, const float sizY, const float sizZ, bool const isSolid) {
  glPushMatrix();
  glTranslatef(begX, begY, begZ);
  glScalef(sizX, sizY, sizZ);
  glTranslatef(0.5f, 0.5f, 0.5f);
  if (isSolid) glutSolidCube(1.0);
  else glutWireCube(1.0);
  glPopMatrix();
}
