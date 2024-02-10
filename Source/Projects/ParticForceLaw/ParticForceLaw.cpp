#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
ParticForceLaw::ParticForceLaw() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void ParticForceLaw::SetActiveProject() {
  if (!isActivProj) {
    D.UI.clear();
    D.UI.push_back(ParamUI("DomainX_________", 1.0));
    D.UI.push_back(ParamUI("DomainY_________", 1.0));
    D.UI.push_back(ParamUI("DomainZ_________", 1.0));
    D.UI.push_back(ParamUI("LatticePitch____", 0.05));
    D.UI.push_back(ParamUI("LatticePattern__", 2));
    D.UI.push_back(ParamUI("ScenarioPreset__", 2));
    D.UI.push_back(ParamUI("ScenarioForce___", 1.0));
    D.UI.push_back(ParamUI("StepsPerDraw____", 4));
    D.UI.push_back(ParamUI("TimeStep________", 0.002));
    D.UI.push_back(ParamUI("IntegType_______", 1));
    D.UI.push_back(ParamUI("ParticleMass____", 1.0));
    D.UI.push_back(ParamUI("VelocityDamping_", 4.0));
    D.UI.push_back(ParamUI("ForceLawPreset__", 0));
    D.UI.push_back(ParamUI("ForceLawNormali_", 1));
    D.UI.push_back(ParamUI("ForceLawScale___", 200.0));
    D.UI.push_back(ParamUI("ForceLawCoeff00_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff01_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff02_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff03_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff04_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff05_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff06_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff07_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff08_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff09_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff10_", 0.0));
    D.UI.push_back(ParamUI("ForceLawCoeff11_", -1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff12_", 0.0));
    D.UI.push_back(ParamUI("ForceLawCoeff13_", 1.0));
    D.UI.push_back(ParamUI("ForceLawCoeff14_", 0.0));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("VisuScale_______", 0.5));
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
bool ParticForceLaw::CheckAlloc() {
  if (D.UI[DomainX_________].hasChanged()) isAllocated= false;
  if (D.UI[DomainY_________].hasChanged()) isAllocated= false;
  if (D.UI[DomainZ_________].hasChanged()) isAllocated= false;
  if (D.UI[LatticePitch____].hasChanged()) isAllocated= false;
  if (D.UI[LatticePattern__].hasChanged()) isAllocated= false;
  if (D.UI[ScenarioPreset__].hasChanged()) isAllocated= false;
  if (D.UI[ScenarioForce___].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool ParticForceLaw::CheckRefresh() {
  if (D.UI[ForceLawPreset__].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawNormali_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawScale___].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff00_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff01_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff02_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff03_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff04_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff05_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff06_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff07_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff08_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff09_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff10_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff11_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff12_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff13_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawCoeff14_].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void ParticForceLaw::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Reset data arrays
  Pos.clear();
  Vel.clear();
  Acc.clear();
  For.clear();
  Col.clear();
  Ext.clear();
  Fix.clear();

  // Get domain dimensions
  D.boxMin= {0.5 - 0.5 * D.UI[DomainX_________].D(), 0.5 - 0.5 * D.UI[DomainY_________].D(), 0.5 - 0.5 * D.UI[DomainZ_________].D()};
  D.boxMax= {0.5 + 0.5 * D.UI[DomainX_________].D(), 0.5 + 0.5 * D.UI[DomainY_________].D(), 0.5 + 0.5 * D.UI[DomainZ_________].D()};

  // Generate the full point cloud over the domain
  std::vector<Vec::Vec3<float>> points;
  float minDist= 0.0f;
  if (D.UI[LatticePattern__].I() == 0) minDist= 1.0f;                    // SCC pattern
  if (D.UI[LatticePattern__].I() == 1) minDist= std::sqrt(3.0f) / 2.0f;  // BCC pattern
  if (D.UI[LatticePattern__].I() == 2) minDist= std::sqrt(2.0f) / 2.0f;  // FCC pattern
  const int nX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * minDist));
  const int nY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * minDist));
  const int nZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * minDist));
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        bool keep= false;
        if (D.UI[LatticePattern__].I() == 0) keep= true;                                                                                // SCC pattern
        if (D.UI[LatticePattern__].I() == 1 && ((x % 2 == 0 && y % 2 == 0 && z % 2 == 0) || (x % 2 + y % 2 + z % 2 == 3))) keep= true;  // BCC pattern
        if (D.UI[LatticePattern__].I() == 2 && ((x + y + z) % 2 == 0)) keep= true;                                                      // FCC pattern
        if (keep)
          points.push_back(Vec::Vec3<float>(
              D.boxMin[0] + x * D.UI[LatticePitch____].F() * minDist,
              D.boxMin[1] + y * D.UI[LatticePitch____].F() * minDist,
              D.boxMin[2] + z * D.UI[LatticePitch____].F() * minDist));
      }
    }
  }

  // Calculate some helper variables for scenario setup
  const Vec::Vec3<float> BoxMin(D.boxMin[0], D.boxMin[1], D.boxMin[2]);
  const Vec::Vec3<float> BoxMax(D.boxMax[0], D.boxMax[1], D.boxMax[2]);
  const float BoxDiag= (BoxMax - BoxMin).norm();

  // Add the subset of points for the current scenario
  for (int k= 0; k < (int)points.size(); k++) {
    const Vec::Vec3<float> RelPos= (points[k] - BoxMin).cwiseDiv(BoxMax - BoxMin);
    // Full set of points
    if (D.UI[ScenarioPreset__].I() == 0) {
      Pos.push_back(points[k]);
      Fix.push_back(false);
    }
    // Box falling on steps
    else if (D.UI[ScenarioPreset__].I() == 1) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.8f)).abs().maxCoeff() < 0.1f * BoxDiag) {
        Pos.push_back(points[k]);
        Ext.push_back(Vec::Vec3<float>(0.0f, 0.0f, -D.UI[ScenarioForce___].F() * D.UI[ParticleMass____].F()));
      }
      else if (RelPos[0] > 0.5f && RelPos[2] > 0.3f && RelPos[2] < 0.5f && (RelPos[1] + RelPos[2]) <= 0.8f && RelPos[1] < 0.5f) {
        Pos.push_back(points[k]);
        Col.push_back(Vec::Vec3<float>(0.3f, 0.3f, 0.3f));
        Fix.push_back(true);
      }
      else if (RelPos[0] < 0.6f && RelPos[2] < 0.06f && RelPos[1] > 0.5f) {
        Pos.push_back(points[k]);
        Col.push_back(Vec::Vec3<float>(0.3f, 0.3f, 0.3f));
        Fix.push_back(true);
      }
    }
    // Ball blasting through wall
    else if (D.UI[ScenarioPreset__].I() == 2) {
      if ((RelPos - Vec::Vec3<float>(0.7f, 0.35f, 0.7f)).norm() < 0.15f) {
        Pos.push_back(points[k]);
        Vel.push_back(Vec::Vec3<float>(0.0f, D.UI[ScenarioForce___].F(), 0.0f));
        Col.push_back(Vec::Vec3<float>(0.8f, 0.3f, 0.3f));
      }
      else if ((RelPos - Vec::Vec3<float>(0.3f, 0.15f, 0.3f)).norm() < 0.15f) {
        Pos.push_back(points[k]);
        Vel.push_back(Vec::Vec3<float>(0.0f, D.UI[ScenarioForce___].F(), 0.0f));
        Col.push_back(Vec::Vec3<float>(0.3f, 0.8f, 0.3f));
      }
      else if (RelPos[1] > 0.7f && RelPos[1] < 0.75f) {
        Pos.push_back(points[k]);
      }
    }
    // Coupon stretch
    else if (D.UI[ScenarioPreset__].I() == 2) {
      if ((RelPos - Vec::Vec3<float>(0.7f, 0.35f, 0.7f)).norm() < 0.15f) {
        Pos.push_back(points[k]);
        Vel.push_back(Vec::Vec3<float>(0.0f, D.UI[ScenarioForce___].F(), 0.0f));
        Col.push_back(Vec::Vec3<float>(0.8f, 0.3f, 0.3f));
      }
      else if ((RelPos - Vec::Vec3<float>(0.3f, 0.15f, 0.3f)).norm() < 0.15f) {
        Pos.push_back(points[k]);
        Vel.push_back(Vec::Vec3<float>(0.0f, D.UI[ScenarioForce___].F(), 0.0f));
        Col.push_back(Vec::Vec3<float>(0.3f, 0.8f, 0.3f));
      }
      else if (RelPos[1] > 0.7f && RelPos[1] < 0.75f) {
        Pos.push_back(points[k]);
      }
    }

    if (Vel.size() < Pos.size()) Vel.resize(Pos.size(), Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Acc.size() < Pos.size()) Acc.resize(Pos.size(), Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (For.size() < Pos.size()) For.resize(Pos.size(), Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Col.size() < Pos.size()) Col.resize(Pos.size(), Vec::Vec3<float>(0.5f, 0.5f, 0.5f));
    if (Ext.size() < Pos.size()) Ext.resize(Pos.size(), Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Fix.size() < Pos.size()) Fix.resize(Pos.size(), false);
  }
}


// Refresh the project
void ParticForceLaw::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Generate the force law
  ForceLaw.clear();
  // Custom force law
  if (D.UI[ForceLawPreset__].I() == 0) {
    ForceLaw.push_back(D.UI[ForceLawCoeff00_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff01_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff02_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff03_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff04_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff05_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff06_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff07_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff08_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff09_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff10_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff11_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff12_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff13_].F());
    ForceLaw.push_back(D.UI[ForceLawCoeff14_].F());
  }
  // Sample force law
  else if (D.UI[ForceLawPreset__].I() == 1) {
    ForceLaw= std::vector<float>{1.5f, 1.5f, 1.5f, 1.5f, 1.5f, 1.5f, 1.5f, 1.5f, 1.5f, 1.5f, 0.0f, -2.2f, 0.0f, 1.2f, 0.0f};
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]*= 1.e8;
  }
  // Elastic force law
  else if (D.UI[ForceLawPreset__].I() == 2) {
    ForceLaw= std::vector<float>{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f};
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]*= 1.e6;
  }
  // Brittle force law
  else if (D.UI[ForceLawPreset__].I() == 3) {
    ForceLaw= std::vector<float>{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, -0.1f, 0.1f, 0.1f, 0.0f};
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]*= 1.e6;
  }
  // Viscous force law
  else if (D.UI[ForceLawPreset__].I() == 4) {
    ForceLaw= std::vector<float>{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.9f, 0.5f, 0.0f, -0.3f, -0.3f, -0.2f, 0.0f};
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]*= 1.e2;
  }
  // Steel AISI 4340 force law
  else if (D.UI[ForceLawPreset__].I() == 5) {
    ForceLaw= std::vector<float>{3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 0.0f, -2.5f, 0.0f, 2.5f, 0.0f};
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]*= 1.e10;
  }

  // Optionally normalize the force law
  if (D.UI[ForceLawNormali_].B()) {
    float refVal= ForceLaw[0];
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]/= refVal;
  }

  // Scale the force law
  for (int k= 0; k < (int)ForceLaw.size(); k++)
    ForceLaw[k]*= D.UI[ForceLawScale___].F();

  // Draw the force law in the scatter plot
  D.Scatter.clear();
  D.Scatter.resize(1);
  D.Scatter[0].name= "ForceLaw";
  for (int k= 0; k < (int)ForceLaw.size(); k++)
    D.Scatter[0].val.push_back(std::array<double, 2>{(double)k / 10.0, (double)ForceLaw[k]});
}


// Handle keypress
void ParticForceLaw::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Animate the project
void ParticForceLaw::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Get and check particle system sizes
  const int nbParticles= (int)Pos.size();
  if (nbParticles <= 0) return;

  // Step forward simulation with explicit numerical integration
  for (int stepIdx= 0; stepIdx < D.UI[StepsPerDraw____].I(); stepIdx++) {
    const float dt= D.UI[TimeStep________].F();
    if (D.UI[IntegType_______].I() == 0) {
      // Evaluate net forces acting on particles
      ComputeForces();
      // Euler integration
      for (int k= 0; k < nbParticles; k++) {
        if (!Fix[k]) {
          Acc[k]= For[k] / D.UI[ParticleMass____].F();  // at+1 = ft / m
          Vel[k]+= Acc[k] * dt;                         // vt+1 = at + at+1 * dt
          Pos[k]+= Vel[k] * dt;                         // xt+1 = vt + vt+1 * dt
        }
      }
    }
    else {
      // Velocity Verlet integration - position update
      for (int k= 0; k < nbParticles; k++) {
        if (!Fix[k]) {
          Pos[k]+= Vel[k] * dt + 0.5 * Acc[k] * dt * dt;  // xt+1 = xt + vt * dt + 0.5 * at * dt * dt
        }
      }
      // Evaluate net forces acting on particles
      ComputeForces();
      // Velocity Verlet integration - acceleration and velocity update
      for (int k= 0; k < nbParticles; k++) {
        if (!Fix[k]) {
          const Vec::Vec3<float> oldAcc= Acc[k];
          Acc[k]= For[k] / D.UI[ParticleMass____].F();  // at+1 = ft+1 / m
          Vel[k]+= 0.5 * (oldAcc + Acc[k]) * dt;        // vt+1 = vt + 0.5 * (at + at+1) * dt
        }
      }
    }
  }

  // Plot data
  D.Plot.resize(1);
  D.Plot[0].name= "KE";
  if (D.Plot[0].val.size() < 10000) {
    D.Plot[0].val.reserve(10000);
    D.Plot[0].val.push_back(0.0f);
    for (int k= 0; k < (int)Pos.size(); k++)
      D.Plot[0].val[D.Plot[0].val.size() - 1]+= D.UI[ParticleMass____].F() * Vel[k].normSquared();
  }
}


// Draw the project
void ParticForceLaw::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Display particles
  if (D.displayMode1) {
    glEnable(GL_LIGHTING);
    for (int k= 0; k < (int)Pos.size(); k++) {
      // Set particle color
      float r= 0.5, g= 0.5, b= 0.5;
      if (D.UI[ColorMode_______].I() == 0) {
        r= Col[k][0];
        g= Col[k][1];
        b= Col[k][2];
      }
      if (D.UI[ColorMode_______].I() == 1) {
        r= D.UI[ColorFactor_____].F() * Ext[k][0] + 0.5f;
        g= D.UI[ColorFactor_____].F() * Ext[k][1] + 0.5f;
        b= D.UI[ColorFactor_____].F() * Ext[k][2] + 0.5f;
      }
      if (D.UI[ColorMode_______].I() == 2) {
        r= D.UI[ColorFactor_____].F() * Vel[k][0] + 0.5f;
        g= D.UI[ColorFactor_____].F() * Vel[k][1] + 0.5f;
        b= D.UI[ColorFactor_____].F() * Vel[k][2] + 0.5f;
      }
      if (D.UI[ColorMode_______].I() == 3) Colormap::RatioToJetBrightSmooth(Vel[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      if (D.UI[ColorMode_______].I() == 4) Colormap::RatioToJetBrightSmooth(For[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      glColor3f(r, g, b);
      // Draw sphere
      glPushMatrix();
      glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
      const float scale= D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F();
      glScalef(scale, scale, scale);
      glutSolidSphere(1.0, 12, 6);
      glPopMatrix();
    }
    glDisable(GL_LIGHTING);
  }
}


void ParticForceLaw::ComputeForces() {
  // Get and check particle system sizes
  const int nbParticles= (int)Pos.size();
  if (nbParticles <= 0) return;

  // Reset forces
  std::fill(For.begin(), For.end(), Vec::Vec3<float>(0.0, 0.0, 0.0));

  // Precompute values
  const float forceReach= D.UI[LatticePitch____].F() * 1.4f;

  // Compute particle forces
  for (int k0= 0; k0 < nbParticles; k0++) {
    // Interaction forces
    for (int k1= k0 + 1; k1 < nbParticles; k1++) {
      // Reject if out of reach
      const Vec::Vec3<float> distVec= Pos[k0] - Pos[k1];
      if (distVec.normSquared() > forceReach * forceReach) continue;
      // Get linear interpolation of force law for the given particle distance
      const float distVal= distVec.norm();
      const float valFloat= (float)(ForceLaw.size() - 1) * distVal / forceReach;
      const int low= std::min(std::max((int)std::floor(valFloat), 0), (int)ForceLaw.size() - 1);
      const int upp= std::min(std::max(low + 1, 0), (int)ForceLaw.size() - 1);
      const float ratio= valFloat - (float)low;
      const float forceVal= (1.0 - ratio) * ForceLaw[low] + (ratio)*ForceLaw[upp];
      // Apply inter-particle force to both particles
      For[k0]+= forceVal * distVec / distVal;
      For[k1]-= forceVal * distVec / distVal;
      // Apply inter-particle damping linearly proportional to difference in radial velocity of particle pair
      const float RadialVel= (Vel[k0] - Vel[k1]).dot(distVec / distVal);
      For[k0]-= D.UI[VelocityDamping_].F() * RadialVel * distVec / distVal;
      For[k1]+= D.UI[VelocityDamping_].F() * RadialVel * distVec / distVal;
    }
    // External forces
    For[k0]+= Ext[k0];
  }
}
