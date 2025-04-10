#include "PosiBasedDynam.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void PosiBasedDynam::DrawScene() {

  // Draw the nodes
  if (D.displayMode[1]) {
    glPointSize(1000.0f * D.UI[RadParticl______].F());
    if (D.UI[VisuSimple______].I() > 0) glBegin(GL_POINTS);
    else                                glEnable(GL_LIGHTING);
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (D.UI[VisuSimple______].I() < 0 && k >= D.UI[NumParticl______].I()) break;
      float r= 0.75f, g= 0.75f, b= 0.75f;
      if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToJetBrightSmooth(Mass[k] * D.UI[ColorFactor_____].F(), r, g, b);
      if (D.UI[ColorMode_______].I() == 2) Colormap::RatioToJetBrightSmooth(Vel[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      glColor3f(r, g, b);
      if (D.UI[VisuSimple______].I() > 0) {
        glVertex3fv(Pos[k].array());
      }
      else {
        glPushMatrix();
        glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
        glutSolidSphere(D.UI[RadParticl______].F(), 12, 6);
        glPopMatrix();
      }

    }
    glPointSize(1.0f);
    if (D.UI[VisuSimple______].I() > 0) glEnd();
    else                                glDisable(GL_LIGHTING);
  }

  // Draw the edges
  if (D.displayMode[2]) {
    glLineWidth(2.0f);
    glColor3f(0.75f, 0.75f, 0.75f);
    glBegin(GL_LINES);
    for (std::array<int, 2> curEdge : Edge) {
      glVertex3fv(Pos[curEdge[0]].array());
      glVertex3fv(Pos[curEdge[1]].array());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the triangles
  if (D.displayMode[3]) {
    glColor3f(0.75f, 0.75f, 0.75f);
    glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (std::array<int, 3> tri : Tri) {
      Vec::Vec3<float> v0(Pos[tri[0]]), v1(Pos[tri[1]]), v2(Pos[tri[2]]);
      Vec::Vec3<float> normal= (v1 - v0).cross(v2 - v0).normalized();
      glNormal3fv(normal.array());
      glVertex3fv(v0.array());
      glVertex3fv(v1.array());
      glVertex3fv(v2.array());
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }

  // Draw the tetrahedras
  if (D.displayMode[4]) {
    glColor3f(0.75f, 0.75f, 0.75f);
    glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (std::array<int, 4> tet : Tet) {
      Vec::Vec3<float> v0, v1, v2;
      const Vec::Vec3<float> center= 0.25f * (Pos[tet[0]] + Pos[tet[1]] + Pos[tet[2]] + Pos[tet[3]]);
      for (int k= 0; k < 4; k++) {
        if      (k == 0) { v0= Pos[tet[0]]; v1= Pos[tet[2]]; v2= Pos[tet[1]]; }
        else if (k == 1) { v0= Pos[tet[0]]; v1= Pos[tet[1]]; v2= Pos[tet[3]]; }
        else if (k == 2) { v0= Pos[tet[0]]; v1= Pos[tet[3]]; v2= Pos[tet[2]]; }
        else if (k == 3) { v0= Pos[tet[1]]; v1= Pos[tet[2]]; v2= Pos[tet[3]]; }
        Vec::Vec3<float> normal= (v1 - v0).cross(v2 - v0).normalized();
        v0= center + (v0 - center) * 0.75f;
        v1= center + (v1 - center) * 0.75f;
        v2= center + (v2 - center) * 0.75f;
        glNormal3fv(normal.array());
        glVertex3fv(v0.array());
        glVertex3fv(v1.array());
        glVertex3fv(v2.array());
      }
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }

  // Draw the triangles normals
  if (D.displayMode[5]) {
    glColor3f(0.75f, 0.75f, 0.75f);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    for (std::array<int, 3> tri : Tri) {
      Vec::Vec3<float> v0(Pos[tri[0]]), v1(Pos[tri[1]]), v2(Pos[tri[2]]);
      Vec::Vec3<float> normal= (v1 - v0).cross(v2 - v0);
      glVertex3fv(((v0 + v1 + v2) / 3.0f).array());
      glVertex3fv(((v0 + v1 + v2) / 3.0f + normal).array());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the mouse ball
  if (D.displayMode[6]) {
    if (D.mouseMiddleButtonState > 0) {
      glColor3f(0.75f, 0.75f, 0.75f);
      glEnable(GL_LIGHTING);
      glPushMatrix();
      glTranslatef(D.mouseProjX[0], D.mouseProjX[1], D.mouseProjX[2]);
      glutSolidSphere(D.UI[RadMouse________].F(), 24, 12);
      glPopMatrix();
      glDisable(GL_LIGHTING);
    }
  }
}
