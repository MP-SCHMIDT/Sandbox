#include "PosiBasedDynam.hpp"


// Standard lib
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
PosiBasedDynam::PosiBasedDynam() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void PosiBasedDynam::SetActiveProject() {
  if (!isActivProj) {
    D.UI.clear();
    D.UI.push_back(ParamUI("NumParticl______", 1000));
    D.UI.push_back(ParamUI("RadParticl______", 0.02));
    D.UI.push_back(ParamUI("DomainX_________", 0.5));
    D.UI.push_back(ParamUI("DomainY_________", 0.5));
    D.UI.push_back(ParamUI("DomainZ_________", 0.3));
    D.UI.push_back(ParamUI("TimeStep________", 0.02));
    D.UI.push_back(ParamUI("VelDecay________", 0.1));
    D.UI.push_back(ParamUI("FactorCondu_____", 2.0));
    D.UI.push_back(ParamUI("ForceGrav_______", -1.0));
    D.UI.push_back(ParamUI("ForceBuoy_______", 2.0));
    D.UI.push_back(ParamUI("HeatInput_______", 0.2));
    D.UI.push_back(ParamUI("HeatOutput______", 0.1));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));
  }

  if (D.UI.size() != VerboseLevel____ + 1) {
    printf("[ERROR] Invalid parameter count in UI\n");
  }

  D.boxMin= {0.0, 0.0, 0.0};
  D.boxMax= {1.0, 1.0, 1.0};

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool PosiBasedDynam::CheckAlloc() {
  if (D.UI[NumParticl______].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool PosiBasedDynam::CheckRefresh() {
  if (D.UI[RadParticl______].hasChanged()) isRefreshed= false;
  if (D.UI[DomainX_________].hasChanged()) isRefreshed= false;
  if (D.UI[DomainY_________].hasChanged()) isRefreshed= false;
  if (D.UI[DomainZ_________].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void PosiBasedDynam::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  N= std::max(D.UI[NumParticl______].I(), 1);

  // Allocate data
  PosOld= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  PosCur= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  VelCur= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  AccCur= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  ForCur= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  ColCur= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  RadCur= std::vector<float>(N, 0.0f);
  MasCur= std::vector<float>(N, 0.0f);
  HotCur= std::vector<float>(N, 0.0f);
}


// Refresh the project
void PosiBasedDynam::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Get domain dimensions
  D.boxMin= {0.5f - D.UI[DomainX_________].F(), 0.5f - D.UI[DomainY_________].F(), 0.5f - D.UI[DomainZ_________].F()};
  D.boxMax= {0.5f + D.UI[DomainX_________].F(), 0.5f + D.UI[DomainY_________].F(), 0.5f + D.UI[DomainZ_________].F()};

  // Initialize with random particle properties
  for (int k= 0; k < N; k++) {
    for (int dim= 0; dim < 3; dim++) {
      PosCur[k][dim]= Random::Val((float)D.boxMin[dim], (float)D.boxMax[dim]);
      ColCur[k][dim]= Random::Val(0.0f, 1.0f);
    }
    RadCur[k]= D.UI[RadParticl______].F();
    MasCur[k]= 1.0f;
    HotCur[k]= Random::Val(0.0f, 1.0f);
  }
  PosOld= PosCur;
}


// Handle keypress
void PosiBasedDynam::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Animate the project
void PosiBasedDynam::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Get UI parameters
  const float dt= D.UI[TimeStep________].F();
  const float velocityDecay= 1.0f - D.UI[VelDecay________].F() * dt;
  const Vec::Vec3<float> vecGrav(0.0f, 0.0f, D.UI[ForceGrav_______].F());
  const Vec::Vec3<float> vecBuoy(0.0f, 0.0f, D.UI[ForceBuoy_______].F());
  const float conductionFactor= D.UI[FactorCondu_____].F();
  const float heatAdd= D.UI[HeatInput_______].F();
  const float heatRem= D.UI[HeatOutput______].F();

  // Add or remove heat to particles based on position in the domain
  for (int k0= 0; k0 < N; k0++) {
    const Vec::Vec3<float> posSource(0.5f * (D.boxMin[0] + D.boxMax[0]), 0.5f * (D.boxMin[1] + D.boxMax[1]), D.boxMin[2]);
    const float radSource= 0.1f * ((D.boxMax[0] - D.boxMin[0]) + (D.boxMax[1] - D.boxMin[1]) + (D.boxMax[2] - D.boxMin[2]));
    if ((posSource - PosCur[k0]).normSquared() < radSource)
      HotCur[k0]+= heatAdd * dt;
    HotCur[k0]-= heatRem * dt;
    HotCur[k0]= std::min(std::max(HotCur[k0], 0.0f), 1.0f);
  }

  // Transfer heat between particles (Gauss Seidel)
  std::vector<float> HotOld= HotCur;
  for (int k0= 0; k0 < N; k0++) {
    for (int k1= k0 + 1; k1 < N; k1++) {
      if ((PosCur[k1] - PosCur[k0]).normSquared() <= 1.1f * (RadCur[k0] + RadCur[k1]) * (RadCur[k0] + RadCur[k1])) {
        float val= conductionFactor * (HotOld[k1] - HotOld[k0]) * dt;
        HotCur[k0]+= val;
        HotCur[k1]-= val;
      }
    }
    HotCur[k0]= std::min(std::max(HotCur[k0], 0.0f), 1.0f);
  }

  // Reset forces
  for (int k0= 0; k0 < N; k0++)
    ForCur[k0].set(0.0f, 0.0f, 0.0f);

  // Add gravity forces
  for (int k0= 0; k0 < N; k0++)
    ForCur[k0]+= vecGrav * MasCur[k0];

  // Add boyancy forces
  for (int k0= 0; k0 < N; k0++)
    ForCur[k0]+= vecBuoy * HotCur[k0];

  // Apply boundary constraint
  for (int k0= 0; k0 < N; k0++)
    for (int dim= 0; dim < 3; dim++)
      PosCur[k0][dim]= std::min(std::max(PosCur[k0][dim], (float)D.boxMin[dim]), (float)D.boxMax[dim]);

  // Apply collision constraint (Gauss Seidel)
  for (int k0= 0; k0 < N; k0++) {
    for (int k1= k0 + 1; k1 < N; k1++) {
      if ((PosCur[k1] - PosCur[k0]).normSquared() <= (RadCur[k0] + RadCur[k1]) * (RadCur[k0] + RadCur[k1])) {
        Vec::Vec3<float> val= (PosCur[k1] - PosCur[k0]).normalized() * 0.5f * ((RadCur[k0] + RadCur[k1]) - (PosCur[k1] - PosCur[k0]).norm());
        PosCur[k0]-= val;
        PosCur[k1]+= val;
      }
    }
  }

  // Deduce velocities
  for (int k0= 0; k0 < N; k0++)
    VelCur[k0]= (PosCur[k0] - PosOld[k0]) / dt;

  // Apply explicit velocity damping
  for (int k0= 0; k0 < N; k0++)
    VelCur[k0]= VelCur[k0] * velocityDecay;

  // Update positions
  PosOld= PosCur;
  for (int k0= 0; k0 < N; k0++) {
    AccCur[k0]= ForCur[k0] / MasCur[k0];
    PosCur[k0]= PosCur[k0] + VelCur[k0] * dt + AccCur[k0] * dt * dt;
  }
}


// Draw the project
void PosiBasedDynam::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  glEnable(GL_LIGHTING);
  for (int k= 0; k < N; k++) {
    glPushMatrix();
    glTranslatef(PosCur[k][0], PosCur[k][1], PosCur[k][2]);
    glScalef(RadCur[k], RadCur[k], RadCur[k]);
    float r, g, b;
    if (D.displayMode1)
      Colormap::RatioToJetSmooth(VelCur[k].norm(), r, g, b);
    else if (D.displayMode2)
      Colormap::RatioToBlackBody(HotCur[k], r, g, b);
    else {
      r= VelCur[k][0];
      g= VelCur[k][1];
      b= VelCur[k][2];
    }
    glColor3f(r, g, b);
    glutSolidSphere(1.0, 32, 16);
    glPopMatrix();
  }
  glDisable(GL_LIGHTING);
}
