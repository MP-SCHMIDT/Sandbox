#include "SpaceTimeWorld.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/Bresenham.hpp"
#include "Type/Field.hpp"
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
SpaceTimeWorld::SpaceTimeWorld() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void SpaceTimeWorld::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("InputFileID_____", 0));
    D.UI.push_back(ParamUI("WorldNbT________", 16));
    D.UI.push_back(ParamUI("WorldNbX________", 50));
    D.UI.push_back(ParamUI("WorldNbY________", 80));
    D.UI.push_back(ParamUI("WorldNbZ________", 80));
    D.UI.push_back(ParamUI("MassReach_______", 8));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("ScreenNbH_______", 100));
    D.UI.push_back(ParamUI("ScreenNbV_______", 100));
    D.UI.push_back(ParamUI("ScreenNbS_______", 50));
    D.UI.push_back(ParamUI("CursorWorldT____", 8));
    D.UI.push_back(ParamUI("TimePersist_____", 0.8));
    D.UI.push_back(ParamUI("FactorCurv______", 1.0));
    D.UI.push_back(ParamUI("FactorDoppl_____", 1.0));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));

    D.displayModeLabel[1]= "Vox";
    D.displayModeLabel[2]= "Flow";
    D.displayModeLabel[3]= "Screen";
    D.displayModeLabel[4]= "Photon";
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool SpaceTimeWorld::CheckAlloc() {
  if (D.UI[InputFileID_____].hasChanged()) isAllocated= false;
  if (D.UI[WorldNbT________].hasChanged()) isAllocated= false;
  if (D.UI[WorldNbX________].hasChanged()) isAllocated= false;
  if (D.UI[WorldNbY________].hasChanged()) isAllocated= false;
  if (D.UI[WorldNbZ________].hasChanged()) isAllocated= false;
  if (D.UI[MassReach_______].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool SpaceTimeWorld::CheckRefresh() {
  if (D.UI[ScreenNbH_______].hasChanged()) isRefreshed= false;
  if (D.UI[ScreenNbV_______].hasChanged()) isRefreshed= false;
  if (D.UI[ScreenNbS_______].hasChanged()) isRefreshed= false;
  if (D.UI[CursorWorldT____].hasChanged()) isRefreshed= false;
  if (D.UI[TimePersist_____].hasChanged()) isRefreshed= false;
  if (D.UI[FactorCurv______].hasChanged()) isRefreshed= false;
  if (D.UI[FactorDoppl_____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void SpaceTimeWorld::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  InitVoxelWorld();
  ComputeWorldFlow();
}


// Refresh the project
void SpaceTimeWorld::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  ComputeScreen();
}


// Handle UI parameter change
void SpaceTimeWorld::ParamChange() {
}


// Handle keypress
void SpaceTimeWorld::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Handle mouse action
void SpaceTimeWorld::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void SpaceTimeWorld::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
}


// Draw the project
void SpaceTimeWorld::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the solid voxels
  if (D.displayMode[1]) {
    int idxT= std::min(std::max(D.UI[CursorWorldT____].I(), 0), worldNbT - 1);
    glPushMatrix();
    glScalef(1.0f / float(worldNbX), 1.0f / float(worldNbY), 1.0f / float(worldNbZ));
    glTranslatef(0.5f, 0.5f, 0.5f);
    for (int x= 0; x < worldNbX; x++) {
      for (int y= 0; y < worldNbY; y++) {
        for (int z= 0; z < worldNbZ; z++) {
          if (worldSolid[idxT][x][y][z]) {
            glPushMatrix();
            glTranslatef(float(x), float(y), float(z));
            glColor3fv(worldColor[idxT][x][y][z].array());
            glutSolidCube(1.0);
            glPopMatrix();
          }
        }
      }
    }
    glPopMatrix();
  }

  // Draw the space time flow field
  if (D.displayMode[2]) {
    int idxT= std::min(std::max(D.UI[CursorWorldT____].I(), 0), worldNbT - 1);
    glBegin(GL_LINES);
    // int displaySkipsize= std::pow((worldNbX * worldNbY * worldNbZ) / 10000, 1.0 / 3.0);
    // for (int x= displaySkipsize / 2; x < worldNbX; x+= displaySkipsize) {
    //   for (int y= displaySkipsize / 2; y < worldNbY; y+= displaySkipsize) {
    //     for (int z= displaySkipsize / 2; z < worldNbZ; z+= displaySkipsize) {
    for (int x= 0; x < worldNbX; x++) {
      for (int y= 0; y < worldNbY; y++) {
        for (int z= 0; z < worldNbZ; z++) {
          // if (worldSolid[idxT][x][y][z]) continue;
          Vec::Vec3<float> flowVec(worldFlows[idxT][x][y][z][1], worldFlows[idxT][x][y][z][2], worldFlows[idxT][x][y][z][3]);
          float r, g, b;
          Colormap::RatioToJetBrightSmooth(0.5 + worldFlows[idxT][x][y][z][0], r, g, b);
          glColor3f(r, g, b);
          Vec::Vec3<float> pos((float(x) + 0.5f) / float(worldNbX), (float(y) + 0.5f) / float(worldNbY), (float(z) + 0.5f) / float(worldNbZ));
          glVertex3fv(pos.array());
          glVertex3fv((pos + D.UI[FactorCurv______].F() * flowVec / float(screenNbS)).array());
        }
      }
    }
    glEnd();
  }

  // Draw the screen
  if (D.displayMode[3]) {
    glPushMatrix();
    glTranslatef(1.0f, 0.0f, 0.0f);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    for (int h= 0; h < screenNbH; h++) {
      for (int v= 0; v < screenNbV; v++) {
        glColor3fv(screenColor[h][v].array());
        glRectf(float(h) / float(screenNbH), float(v) / float(screenNbV), float(h + 1) / float(screenNbH), float(v + 1) / float(screenNbV));
      }
    }
    glPopMatrix();
  }

  // Draw the photon paths
  if (D.displayMode[4]) {
    glBegin(GL_LINES);
    // for (int h= 0; h < screenNbH; h++) {
    //   for (int v= 0; v < screenNbV; v++) {
    int displaySkipsize= std::sqrt((screenNbH * screenNbV) / 400);
    for (int h= displaySkipsize / 2; h < screenNbH; h+= displaySkipsize) {
      for (int v= displaySkipsize / 2; v < screenNbV; v+= displaySkipsize) {
        for (int s= 0; s < screenCount[h][v] - 1; s++) {
          Vec::Vec3<float> photonBeg(photonPos[h][v][s][1], photonPos[h][v][s][2], photonPos[h][v][s][3]);
          Vec::Vec3<float> photonEnd(photonPos[h][v][s + 1][1], photonPos[h][v][s + 1][2], photonPos[h][v][s + 1][3]);
          glColor3fv(screenColor[h][v].array());
          glVertex3fv(photonBeg.array());
          glColor3fv(screenColor[h][v].array());
          glVertex3fv(photonEnd.array());
        }
      }
    }
    glEnd();
    glPointSize(2.0f);
    glBegin(GL_POINTS);
    // for (int h= 0; h < screenNbH; h++) {
    //   for (int v= 0; v < screenNbV; v++) {
    for (int h= displaySkipsize / 2; h < screenNbH; h+= displaySkipsize) {
      for (int v= displaySkipsize / 2; v < screenNbV; v+= displaySkipsize) {
        for (int s= 0; s < screenCount[h][v]; s++) {
          Vec::Vec3<float> photonBeg(photonPos[h][v][s][1], photonPos[h][v][s][2], photonPos[h][v][s][3]);
          glColor3fv(screenColor[h][v].array());
          glVertex3fv(photonBeg.array());
        }
      }
    }
    glEnd();
    glPointSize(1.0f);
  }
}
