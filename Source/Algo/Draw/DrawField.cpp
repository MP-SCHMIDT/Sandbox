#include "DrawField.hpp"

// Standard lib
#include <array>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Geom/BoxGrid.hpp"
#include "Type/Field.hpp"


// Draw a scalar field as solid voxel boundary or as layered panels for transparency
void DrawField::DrawColored3DField(const Field::Field3<char>& iShow,
                                   const Field::Field3<std::array<float, 4>>& iColor,
                                   const std::array<double, 3>& iBoxMin,
                                   const std::array<double, 3>& iBoxMax,
                                   const std::array<double, 3>& iCamDir,
                                   const bool iWireframe,
                                   const bool iShading,
                                   const bool iAlphaPanels) {
  // Get dimensions from the fields and box
  if (iShow.nX <= 0 || iShow.nZ <= 0 || iShow.nZ <= 0) return;
  if (iShow.nX != iColor.nX || iShow.nY != iColor.nY || iShow.nZ != iColor.nZ) return;
  const int nX= iShow.nX;
  const int nY= iShow.nY;
  const int nZ= iShow.nZ;
  double voxSizeX= 1.0, voxSizeY= 1.0, voxSizeZ= 1.0;
  BoxGrid::GetVoxelSizes(iShow.nX, iShow.nY, iShow.nZ, iBoxMin, iBoxMax, true, voxSizeX, voxSizeY, voxSizeZ);

  // Toggle wireframe or fill mode
  if (iWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Precompute the panels direction and ordering strategy if in alpha transparency mode
  int panelsMode= 0;
  if (iAlphaPanels) {
    if (std::abs(iCamDir[0]) > std::abs(iCamDir[1]) && std::abs(iCamDir[0]) > std::abs(iCamDir[2])) panelsMode= (iCamDir[0] > 0.0) ? 1 : -1;
    if (std::abs(iCamDir[1]) > std::abs(iCamDir[0]) && std::abs(iCamDir[1]) > std::abs(iCamDir[2])) panelsMode= (iCamDir[1] > 0.0) ? 2 : -2;
    if (std::abs(iCamDir[2]) > std::abs(iCamDir[0]) && std::abs(iCamDir[2]) > std::abs(iCamDir[1])) panelsMode= (iCamDir[2] > 0.0) ? 3 : -3;
  }

  // Set the initial transformation
  glPushMatrix();
  glTranslatef(iBoxMin[0] + 0.5f * voxSizeX, iBoxMin[1] + 0.5f * voxSizeY, iBoxMin[2] + 0.5f * voxSizeZ);
  glScalef(voxSizeX, voxSizeY, voxSizeZ);

  // Draw polygons
  if (iShading) glEnable(GL_LIGHTING);
  glBegin(GL_QUADS);
  // Nested loops with half of them inactive to allow arbitrary sweep ordering
  for (int xMajor= 0; xMajor < ((std::abs(panelsMode) == 1) ? nX : 1); xMajor++) {
    for (int yMajor= 0; yMajor < ((std::abs(panelsMode) == 2) ? nY : 1); yMajor++) {
      for (int zMajor= 0; zMajor < ((std::abs(panelsMode) == 3) ? nZ : 1); zMajor++) {
        for (int xMinor= 0; xMinor < ((std::abs(panelsMode) != 1) ? nX : 1); xMinor++) {
          for (int yMinor= 0; yMinor < ((std::abs(panelsMode) != 2) ? nY : 1); yMinor++) {
            for (int zMinor= 0; zMinor < ((std::abs(panelsMode) != 3) ? nZ : 1); zMinor++) {
              // Determine the sweep indices
              const int xInc= (std::abs(panelsMode) == 1) ? xMajor : xMinor;
              const int yInc= (std::abs(panelsMode) == 2) ? yMajor : yMinor;
              const int zInc= (std::abs(panelsMode) == 3) ? zMajor : zMinor;
              const int x= (std::abs(panelsMode) == 1 && iCamDir[0] < 0.0) ? nX - 1 - xInc : xInc;
              const int y= (std::abs(panelsMode) == 2 && iCamDir[1] < 0.0) ? nY - 1 - yInc : yInc;
              const int z= (std::abs(panelsMode) == 3 && iCamDir[2] < 0.0) ? nZ - 1 - zInc : zInc;
              // Draw the polygons based on visibility and panel mode
              if (!iShow.at(x, y, z)) continue;
              glColor4f(iColor.at(x, y, z)[0], iColor.at(x, y, z)[1], iColor.at(x, y, z)[2], iColor.at(x, y, z)[3]);
              if (iAlphaPanels) {
                if (std::abs(panelsMode) == 1) {
                  glNormal3f(1.0f, 0.0f, 0.0f);
                  glVertex3f(x, y - 0.5f, z - 0.5f);
                  glVertex3f(x, y + 0.5f, z - 0.5f);
                  glVertex3f(x, y + 0.5f, z + 0.5f);
                  glVertex3f(x, y - 0.5f, z + 0.5f);
                }
                if (std::abs(panelsMode) == 2) {
                  glNormal3f(0.0f, 1.0f, 0.0f);
                  glVertex3f(x - 0.5f, y, z - 0.5f);
                  glVertex3f(x + 0.5f, y, z - 0.5f);
                  glVertex3f(x + 0.5f, y, z + 0.5f);
                  glVertex3f(x - 0.5f, y, z + 0.5f);
                }
                if (std::abs(panelsMode) == 3) {
                  glNormal3f(0.0f, 0.0f, 1.0f);
                  glVertex3f(x - 0.5f, y - 0.5f, z);
                  glVertex3f(x + 0.5f, y - 0.5f, z);
                  glVertex3f(x + 0.5f, y + 0.5f, z);
                  glVertex3f(x - 0.5f, y + 0.5f, z);
                }
              }
              else {
                if (x == 0 || !iShow.at(x - 1, y, z)) {
                  glNormal3f(-1.0f, 0.0f, 0.0f);
                  glVertex3f(x - 0.5f, y - 0.5f, z - 0.5f);
                  glVertex3f(x - 0.5f, y + 0.5f, z - 0.5f);
                  glVertex3f(x - 0.5f, y + 0.5f, z + 0.5f);
                  glVertex3f(x - 0.5f, y - 0.5f, z + 0.5f);
                }
                if (x == nX - 1 || !iShow.at(x + 1, y, z)) {
                  glNormal3f(+1.0f, 0.0f, 0.0f);
                  glVertex3f(x + 0.5f, y - 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y + 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y + 0.5f, z + 0.5f);
                  glVertex3f(x + 0.5f, y - 0.5f, z + 0.5f);
                }
                if (y == 0 || !iShow.at(x, y - 1, z)) {
                  glNormal3f(0.0f, -1.0f, 0.0f);
                  glVertex3f(x - 0.5f, y - 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y - 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y - 0.5f, z + 0.5f);
                  glVertex3f(x - 0.5f, y - 0.5f, z + 0.5f);
                }
                if (y == nY - 1 || !iShow.at(x, y + 1, z)) {
                  glNormal3f(0.0f, +1.0f, 0.0f);
                  glVertex3f(x - 0.5f, y + 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y + 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y + 0.5f, z + 0.5f);
                  glVertex3f(x - 0.5f, y + 0.5f, z + 0.5f);
                }
                if (z == 0 || !iShow.at(x, y, z - 1)) {
                  glNormal3f(0.0f, 0.0f, -1.0f);
                  glVertex3f(x - 0.5f, y - 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y - 0.5f, z - 0.5f);
                  glVertex3f(x + 0.5f, y + 0.5f, z - 0.5f);
                  glVertex3f(x - 0.5f, y + 0.5f, z - 0.5f);
                }
                if (z == nZ - 1 || !iShow.at(x, y, z + 1)) {
                  glNormal3f(0.0f, 0.0f, +1.0f);
                  glVertex3f(x - 0.5f, y - 0.5f, z + 0.5f);
                  glVertex3f(x + 0.5f, y - 0.5f, z + 0.5f);
                  glVertex3f(x + 0.5f, y + 0.5f, z + 0.5f);
                  glVertex3f(x - 0.5f, y + 0.5f, z + 0.5f);
                }
              }
            }
          }
        }
      }
    }
  }
  glEnd();
  if (iShading) glDisable(GL_LIGHTING);
  glPopMatrix();
}
