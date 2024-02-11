#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

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
    D.UI.push_back(ParamUI("ScenarioPreset__", 5));
    D.UI.push_back(ParamUI("Scenario2DID____", 0));
    D.UI.push_back(ParamUI("Scenario3DID____", 0));
    D.UI.push_back(ParamUI("Scenario2DThick_", 0.1));
    D.UI.push_back(ParamUI("LatticePitch____", 0.05));
    D.UI.push_back(ParamUI("LatticePattern__", 2));
    D.UI.push_back(ParamUI("BCVelX__________", 0.0));
    D.UI.push_back(ParamUI("BCVelY__________", 0.0));
    D.UI.push_back(ParamUI("BCVelZ__________", 1.0));
    D.UI.push_back(ParamUI("BCForX__________", 0.0));
    D.UI.push_back(ParamUI("BCForY__________", 0.0));
    D.UI.push_back(ParamUI("BCForZ__________", 1.0));
    D.UI.push_back(ParamUI("StepsPerDraw____", 1));
    D.UI.push_back(ParamUI("TimeStep________", 0.002));
    D.UI.push_back(ParamUI("IntegType_______", 1));
    D.UI.push_back(ParamUI("ParticleMass____", 1.0));
    D.UI.push_back(ParamUI("VelocityDamping_", 1.0));
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
    D.UI.push_back(ParamUI("ColorMode_______", 2));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("VisuScale_______", 0.5));
    D.UI.push_back(ParamUI("VisuSimple______", 1));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));
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
  if (D.UI[ScenarioPreset__].hasChanged()) isAllocated= false;
  if (D.UI[Scenario2DID____].hasChanged()) isAllocated= false;
  if (D.UI[Scenario3DID____].hasChanged()) isAllocated= false;
  if (D.UI[Scenario2DThick_].hasChanged()) isAllocated= false;
  if (D.UI[LatticePitch____].hasChanged()) isAllocated= false;
  if (D.UI[LatticePattern__].hasChanged()) isAllocated= false;
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
  BCPos.clear();
  BCVel.clear();
  BCFor.clear();

  // Get domain dimensions
  D.boxMin= {0.5 - 0.5 * D.UI[DomainX_________].D(), 0.5 - 0.5 * D.UI[DomainY_________].D(), 0.5 - 0.5 * D.UI[DomainZ_________].D()};
  D.boxMax= {0.5 + 0.5 * D.UI[DomainX_________].D(), 0.5 + 0.5 * D.UI[DomainY_________].D(), 0.5 + 0.5 * D.UI[DomainZ_________].D()};

  // Generate the full point cloud over the domain
  std::vector<Vec::Vec3<float>> pointCloud;
  BuildBaseCloud(pointCloud);

  // Generate the scenario
  BuildScenario(pointCloud);
}


// Refresh the project
void ParticForceLaw::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Generate the force law
  BuildForceLaw();

  // Draw the force law in the plot
  if (D.Plot.size() < 2) D.Plot.resize(2);
  D.Plot[1].name= "ForceLaw";
  D.Plot[1].val.resize(ForceLaw.size());
  for (int k= 0; k < (int)ForceLaw.size(); k++)
    D.Plot[1].val[k]= ForceLaw[k];
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

  // Step forward simulation with explicit numerical integration
  for (int stepIdx= 0; stepIdx < D.UI[StepsPerDraw____].I(); stepIdx++) {
    StepSimulation();
  }

  // Plot data
  if (D.Plot.size() < 1) D.Plot.resize(1);
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

  if (D.UI[VerboseLevel____].I() >= 2) Timer::PushTimer();

  // Display particles
  if (D.displayMode1) {
    if (D.UI[VisuSimple______].B()) {
      glPointSize(1000.0f * D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F());
      glBegin(GL_POINTS);
    }
    else {
      glEnable(GL_LIGHTING);
    }
    for (int k= 0; k < (int)Pos.size(); k++) {
      // Set particle color
      float r= 0.5, g= 0.5, b= 0.5;
      if (D.UI[ColorMode_______].I() == 0) {
        r= Col[k][0];
        g= Col[k][1];
        b= Col[k][2];
      }
      if (D.UI[ColorMode_______].I() == 1) {
        r= 0.5f + 0.3f * (float)BCPos[k];
        g= 0.5f + 0.3f * (float)BCVel[k];
        b= 0.5f + 0.3f * (float)BCFor[k];
      }
      if (D.UI[ColorMode_______].I() == 2) {
        r= D.UI[ColorFactor_____].F() * Vel[k][0] + 0.5f;
        g= D.UI[ColorFactor_____].F() * Vel[k][1] + 0.5f;
        b= D.UI[ColorFactor_____].F() * Vel[k][2] + 0.5f;
      }
      if (D.UI[ColorMode_______].I() == 3) Colormap::RatioToJetBrightSmooth(Vel[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      if (D.UI[ColorMode_______].I() == 4) Colormap::RatioToJetBrightSmooth(For[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      glColor3f(r, g, b);
      if (D.UI[VisuSimple______].B()) {
        glVertex3fv(Pos[k].array());
      }
      else {
        glPushMatrix();
        glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
        const float scale= D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F();
        glScalef(scale, scale, scale);
        glutSolidSphere(1.0, 12, 6);
        glPopMatrix();
      }
    }
    if (D.UI[VisuSimple______].B()) {
      glEnd();
    }
    else {
      glDisable(GL_LIGHTING);
    }
  }
  if (D.UI[VerboseLevel____].I() >= 2) printf("DrawT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  oPointCloud.clear();
  float minDist= 0.0f;
  if (D.UI[LatticePattern__].I() == 0) minDist= 1.0f;                                    // SCC pattern
  if (D.UI[LatticePattern__].I() == 1) minDist= (2.0f / 3.0f) * std::sqrt(3.0f) / 2.0f;  // BCC pattern
  if (D.UI[LatticePattern__].I() == 2) minDist= std::sqrt(2.0f) / 2.0f;                  // FCC pattern
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
          oPointCloud.push_back(Vec::Vec3<float>(
              D.boxMin[0] + x * D.UI[LatticePitch____].F() * minDist,
              D.boxMin[1] + y * D.UI[LatticePitch____].F() * minDist,
              D.boxMin[2] + z * D.UI[LatticePitch____].F() * minDist));
      }
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("BaseCloudT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildScenario(const std::vector<Vec::Vec3<float>>& iPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Calculate some helper variables for scenario setup
  const Vec::Vec3<float> BoxMin(D.boxMin[0], D.boxMin[1], D.boxMin[2]);
  const Vec::Vec3<float> BoxMax(D.boxMax[0], D.boxMax[1], D.boxMax[2]);
  const float BoxDiag= (BoxMax - BoxMin).norm();
  const Vec::Vec3<float> BCForVecNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCForVecPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCVelVecNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> BCVelVecPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());

  // Load the 2D scenario file if needed
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (D.UI[ScenarioPreset__].I() == 0) {
    if (D.UI[Scenario2DID____].I() == 0) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/Coupon.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 1) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/SandiaFracture.bmp", imageRGBA, false);
  }

  // Load the 3D scenario file if needed
  if (D.UI[ScenarioPreset__].I() == 1) {
    // TODO
  }

  // Add the subset of points for the current scenario
  for (int k= 0; k < (int)iPointCloud.size(); k++) {
    const Vec::Vec3<float> RelPos= (iPointCloud[k] - BoxMin).cwiseDiv(BoxMax - BoxMin);
    // 2D loaded scenario
    if (D.UI[ScenarioPreset__].I() == 0) {
      if (!imageRGBA.empty()) {
        if (std::abs(RelPos[0] - 0.5f) < 0.5f * D.UI[Scenario2DThick_].F()) {
          const float posW= (float)(imageRGBA.size() - 1) * RelPos[1];
          const float posH= (float)(imageRGBA[0].size() - 1) * RelPos[2];
          const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
          const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
          const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
          if (colRGBA[3] > 0.5f) {
            Pos.push_back(iPointCloud[k]);
            if (colRGBA[0] > 0.9f) BCPos.push_back(1);
            else if (colRGBA[1] < 0.1f) BCVel.push_back(-1);
            else if (colRGBA[1] > 0.9f) BCVel.push_back(1);
            else if (colRGBA[2] < 0.1f) BCFor.push_back(-1);
            else if (colRGBA[2] > 0.9f) BCFor.push_back(1);
          }
        }
      }
    }
    // 3D loaded scenario
    else if (D.UI[ScenarioPreset__].I() == 1) {
      // TODO set particle based on loaded mesh
    }
    // Full set of points
    else if (D.UI[ScenarioPreset__].I() == 2) {
      Pos.push_back(iPointCloud[k]);
    }
    // Box falling on steps
    else if (D.UI[ScenarioPreset__].I() == 3) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.8f)).abs().maxCoeff() < 0.1f * BoxDiag) {
        Pos.push_back(iPointCloud[k]);
        BCFor.push_back(-1);
      }
      else if (RelPos[0] > 0.5f && RelPos[2] > 0.3f && RelPos[2] < 0.5f && (RelPos[1] + RelPos[2]) <= 0.8f && RelPos[1] < 0.5f) {
        Pos.push_back(iPointCloud[k]);
        Col.push_back(Vec::Vec3<float>(0.3f, 0.3f, 0.3f));
        BCPos.push_back(1);
      }
      else if (RelPos[0] < 0.6f && RelPos[2] < 0.06f && RelPos[1] > 0.5f) {
        Pos.push_back(iPointCloud[k]);
        Col.push_back(Vec::Vec3<float>(0.3f, 0.3f, 0.3f));
        BCPos.push_back(1);
      }
    }
    // Balls blasting through wall
    else if (D.UI[ScenarioPreset__].I() == 4) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.65f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Col.push_back(Vec::Vec3<float>(0.8f, 0.3f, 0.3f));
      }
      else if ((RelPos - Vec::Vec3<float>(0.5f, 0.6f, 0.15f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Col.push_back(Vec::Vec3<float>(0.3f, 0.8f, 0.3f));
      }
      else if (RelPos[2] > 0.9f && RelPos[2] < 0.95f) {
        Pos.push_back(iPointCloud[k]);
      }
    }
    // Coupon stretch - velocity driven
    else if (D.UI[ScenarioPreset__].I() == 5) {
      if (RelPos[0] > 0.45f && RelPos[0] < 0.55f) {
        if (RelPos[1] > 0.3f && RelPos[1] < 0.7f) {
          Pos.push_back(iPointCloud[k]);
          if (RelPos[2] < 0.1f) BCVel.push_back(-1);
          else if (RelPos[2] > 0.9f) BCVel.push_back(1);
        }
      }
    }
    // Coupon stretch - force driven
    else if (D.UI[ScenarioPreset__].I() == 6) {
      if (RelPos[0] > 0.45f && RelPos[0] < 0.55f) {
        if (RelPos[1] > 0.3f && RelPos[1] < 0.7f) {
          Pos.push_back(iPointCloud[k]);
          if (RelPos[2] < 0.1f) BCFor.push_back(-1);
          else if (RelPos[2] > 0.9f) BCFor.push_back(1);
        }
      }
    }

    if (BCPos.size() < Pos.size()) BCPos.push_back(0);
    if (BCVel.size() < Pos.size()) BCVel.push_back(0);
    if (BCFor.size() < Pos.size()) BCFor.push_back(0);

    if (Vel.size() < Pos.size()) Vel.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Acc.size() < Pos.size()) Acc.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (For.size() < Pos.size()) For.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Col.size() < Pos.size()) Col.push_back(Vec::Vec3<float>(0.5f, 0.5f, 0.5f));

    if (!Pos.empty()) {
      if (BCVel[Pos.size() - 1] != 0) Vel[Pos.size() - 1]= (BCVel[Pos.size() - 1] < 0) ? (BCVelVecNega) : (BCVelVecPosi);
      if (BCFor[Pos.size() - 1] != 0) For[Pos.size() - 1]= (BCFor[Pos.size() - 1] < 0) ? (BCForVecNega) : (BCForVecPosi);
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("ScenarioT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildForceLaw() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

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

  if (D.UI[VerboseLevel____].I() >= 1) printf("ForceLawT %f\n", Timer::PopTimer());
}


void ParticForceLaw::ComputeForces() {
  // Get and check particle system sizes
  const int nbParticles= (int)Pos.size();
  if (nbParticles <= 0) return;

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset forces
  std::fill(For.begin(), For.end(), Vec::Vec3<float>(0.0, 0.0, 0.0));

  // Precompute values
  const float forceReach= D.UI[LatticePitch____].F() * 1.4f;
  const Vec::Vec3<float> BCForVecNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCForVecPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());

  // Compute particle forces
  for (int k0= 0; k0 < nbParticles; k0++) {
    // Interaction forces
    // TODO spatial partition buckets
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
    if (BCFor[k0] != 0) For[k0]+= (BCFor[k0] < 0) ? (BCForVecNega) : (BCForVecPosi);
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("ForcesT %f\n", Timer::PopTimer());
}


void ParticForceLaw::StepSimulation() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Precompute values
  const Vec::Vec3<float> BCVelVecNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> BCVelVecPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());

  const float dt= D.UI[TimeStep________].F();
  if (D.UI[IntegType_______].I() == 0) {
    // Evaluate net forces acting on particles
    ComputeForces();
    // Euler integration
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (BCPos[k] == 0) {
        Acc[k]= For[k] / D.UI[ParticleMass____].F();                // at+1 = ft / m
        if (BCVel[k] == 0) Vel[k]+= Acc[k] * dt;                    // vt+1 = at + at+1 * dt
        else Vel[k]= (BCVel[k] < 0) ? BCVelVecNega : BCVelVecPosi;  // Prescribed velocity
        Pos[k]+= Vel[k] * dt;                                       // xt+1 = vt + vt+1 * dt
      }
    }
  }
  else {
    // Velocity Verlet integration - position update
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (BCPos[k] == 0) {
        Pos[k]+= Vel[k] * dt + 0.5f * Acc[k] * dt * dt;  // xt+1 = xt + vt * dt + 0.5 * at * dt * dt
      }
    }
    // Evaluate net forces acting on particles
    ComputeForces();
    // Velocity Verlet integration - acceleration and velocity update
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (BCPos[k] == 0) {
        const Vec::Vec3<float> oldAcc= Acc[k];
        Acc[k]= For[k] / D.UI[ParticleMass____].F();                // at+1 = ft+1 / m
        if (BCVel[k] == 0) Vel[k]+= 0.5f * (oldAcc + Acc[k]) * dt;  // vt+1 = vt + 0.5 * (at + at+1) * dt
        else Vel[k]= (BCVel[k] < 0) ? BCVelVecNega : BCVelVecPosi;  // Prescribed velocity
      }
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("StepT %f\n", Timer::PopTimer());
}
