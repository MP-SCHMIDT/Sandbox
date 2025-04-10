#include "NBodyGravDynam.hpp"


// Standard lib
#include <cmath>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void NBodyGravDynam::DrawScene() {
  // Draw the bodies
  if (D.displayMode[1]) {
    glPointSize(1000.0f * D.UI[BodyRadius______].F());
    if (D.UI[SphereSimple____].I() == 1) glEnable(GL_POINT_SMOOTH);
    if (D.UI[SphereSimple____].I() > 0)  glBegin(GL_POINTS);
    else                                 glEnable(GL_LIGHTING);
    for (unsigned int k0= 0; k0 < N; k0++) {
      float r, g, b;
      if (D.UI[ColorMode_______].I() == 0)      Colormap::RatioToJetBrightSmooth(Vel[k0].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      else if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToPlasma(std::sqrt(For[k0].norm()) * D.UI[ColorFactor_____].F(), r, g, b);
      else if (D.UI[ColorMode_______].I() == 2) Colormap::RatioToRainbow(float(k0)/float(N), r, g, b);
      glColor3f(r, g, b);

      if (D.UI[SphereSimple____].I() > 0) {
        glVertex3fv(Pos[k0].array());
      }
      else {
        glPushMatrix();
        glTranslatef(Pos[k0][0], Pos[k0][1], Pos[k0][2]);
        glutSolidSphere(D.UI[BodyRadius______].F(), 12, 6);
        glPopMatrix();
      }
    }
    if (D.UI[SphereSimple____].I() > 0)  glEnd();
    else                                 glDisable(GL_LIGHTING);
    if (D.UI[SphereSimple____].I() == 1) glDisable(GL_POINT_SMOOTH);
    glPointSize(1.0f);
  }
  
  // Draw the bodies velocities
  if (D.displayMode[2]) {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    for (unsigned int k0= 0; k0 < N; k0++) {
      const float r= 0.5f + Vel[k0][0] * D.UI[ColorFactor_____].F();
      const float g= 0.5f + Vel[k0][1] * D.UI[ColorFactor_____].F();
      const float b= 0.5f + Vel[k0][2] * D.UI[ColorFactor_____].F();
      glColor3f(r, g, b);
      glVertex3fv(Pos[k0].array());
      glVertex3fv((Pos[k0] + Vel[k0] * D.UI[ScaleFactor_____].F()).array());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the traversal order through the body list
  if (D.displayMode[3]) {
    glLineWidth(2.0f);
    glBegin(GL_LINE_STRIP);
    for (unsigned int k0= 0; k0 < N; k0++) {
      float r, g, b;
      Colormap::RatioToRainbow(float(k0)/float(N), r, g, b);
      glColor3f(r, g, b);
      glVertex3fv(Pos[k0].array());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the octree
  if (D.displayMode[4]) {
    for (NBodyGravDynam::OctreeNode Cell : Tree) {
      if (D.UI[ShowEmptyCells__].I() || Cell.Count > 0) {
        float r, g, b;
        Colormap::RatioToRainbow(0.1f * Tree[0].Size / Cell.Size, r, g, b);
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(Cell.Center[0], Cell.Center[1], Cell.Center[2]);
        glutWireCube(Cell.Size);
        glPopMatrix();
      }
    }
  }

  // Draw the octree average positions
  if (D.displayMode[5]) {
    glLineWidth(2.0f);
    for (NBodyGravDynam::OctreeNode Cell : Tree) {
      if (Cell.Count <= 0) continue;
      float r, g, b;
      Colormap::RatioToTurbo((float)Cell.Count * D.UI[ColorFactor_____].F(), r, g, b);
      glColor3f(r, g, b);
      glPushMatrix();
      glTranslatef(Cell.AvgPos[0], Cell.AvgPos[1], Cell.AvgPos[2]);
      glScalef(0.1f * Cell.Size, 0.1f * Cell.Size, 0.1f * Cell.Size);
      glutWireIcosahedron();
      glPopMatrix();
    }
    glLineWidth(1.0f);
  }

  // Draw the octree average velocities
  if (D.displayMode[6]) {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    for (NBodyGravDynam::OctreeNode Cell : Tree) {
      if (Cell.Count <= 0) continue;
      const float r= 0.5f + Cell.AvgVel[0] * D.UI[ColorFactor_____].F();
      const float g= 0.5f + Cell.AvgVel[1] * D.UI[ColorFactor_____].F();
      const float b= 0.5f + Cell.AvgVel[2] * D.UI[ColorFactor_____].F();
      glColor3f(r, g, b);
      glVertex3fv(Cell.AvgPos.array());
      glVertex3fv((Cell.AvgPos + Cell.AvgVel * D.UI[ScaleFactor_____].F()).array());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the octree traversal order through the array
  if (D.displayMode[7]) {
    glLineWidth(2.0f);
    glBegin(GL_LINE_STRIP);
    for (unsigned int idxCell= 0; idxCell < Tree.size(); idxCell++) {
      if (Tree[idxCell].Count == 0) continue;
      float r, g, b;
      Colormap::RatioToRainbow(float(idxCell)/float(Tree.size()), r, g, b);
      glColor3f(r, g, b);
      glVertex3fv(Tree[idxCell].Center.array());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the contributing forces for the chosen particle
  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  if (D.displayMode[8]) {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    for (unsigned int k= 0; k < ContribPos.size(); k++) {
      float r, g, b;
      Colormap::RatioToTurbo((float)ContribCount[k] * D.UI[ColorFactor_____].F(), r, g, b);
      glColor3f(r, g, b);
      glVertex3fv(Pos[N-1].array());
      glVertex3fv(ContribPos[k].array());
    }
    glEnd();
    glLineWidth(1.0f);
  }
  #endif
  
  // Draw the contributing cells for the chosen particle
  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  if (D.displayMode[9]) {
    for (unsigned int idxCell : ContribCell) {
      const NBodyGravDynam::OctreeNode Cell= Tree[idxCell];
      float r, g, b;
      Colormap::RatioToRainbow(0.1f * Tree[0].Size / Cell.Size, r, g, b);
      glColor3f(r, g, b);
      glPushMatrix();
      glTranslatef(Cell.Center[0], Cell.Center[1], Cell.Center[2]);
      glutWireCube(Cell.Size);
      glPopMatrix();
    }
  }
  #endif
}
