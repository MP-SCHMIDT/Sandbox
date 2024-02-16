#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/Sketch.hpp"
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
    D.UI.push_back(ParamUI("ScenarioPreset__", 0));
    D.UI.push_back(ParamUI("Scenario2DID____", 0));
    D.UI.push_back(ParamUI("Scenario2DThick_", 0.5));
    D.UI.push_back(ParamUI("LatticePitch____", 0.04));
    D.UI.push_back(ParamUI("MaterialDensity_", 10000.0));
    D.UI.push_back(ParamUI("SpatialSortNb___", 1000));
    D.UI.push_back(ParamUI("SpatialSortSize_", 40));
    D.UI.push_back(ParamUI("BCVelX__________", 0.0));
    D.UI.push_back(ParamUI("BCVelY__________", 0.0));
    D.UI.push_back(ParamUI("BCVelZ__________", 1.0));
    D.UI.push_back(ParamUI("BCForX__________", 0.0));
    D.UI.push_back(ParamUI("BCForY__________", 0.0));
    D.UI.push_back(ParamUI("BCForZ__________", 1.0));
    D.UI.push_back(ParamUI("StepsPerDraw____", 1));
    D.UI.push_back(ParamUI("TimeStep________", 0.002));
    D.UI.push_back(ParamUI("IntegType_______", 1));
    D.UI.push_back(ParamUI("UseForceControl_", 0));
    D.UI.push_back(ParamUI("BCPosCoeff______", 1.0));
    D.UI.push_back(ParamUI("BCVelCoeff______", 1.0));
    D.UI.push_back(ParamUI("VelocityDamping_", 0.0));
    D.UI.push_back(ParamUI("ColorMode_______", 3));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("ColorDecay______", 0.5));
    D.UI.push_back(ParamUI("VisuScale_______", 0.5));
    D.UI.push_back(ParamUI("VisuSimple______", 1));
    D.UI.push_back(ParamUI("ForceLawPreset__", 0));
    D.UI.push_back(ParamUI("ForceLawTailTol_", 1.e-2));
    D.UI.push_back(ParamUI("ForceLawNormali_", 1));
    D.UI.push_back(ParamUI("ForceLawScale___", 100.0));
    D.UI.push_back(ParamUI("ForceLawA_______", 1.0));
    D.UI.push_back(ParamUI("ForceLaw08______", 1.0));
    D.UI.push_back(ParamUI("ForceLaw09______", 1.0));
    D.UI.push_back(ParamUI("ForceLawB_______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw11______", -1.32));
    D.UI.push_back(ParamUI("ForceLaw12______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw13______", 1.32));
    D.UI.push_back(ParamUI("ForceLawC_______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw15______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw20______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw25______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw30______", 0.0));
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
  if (D.UI[ScenarioPreset__].hasChanged()) isAllocated= false;
  if (D.UI[Scenario2DID____].hasChanged()) isAllocated= false;
  if (D.UI[Scenario2DThick_].hasChanged()) isAllocated= false;
  if (D.UI[LatticePitch____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool ParticForceLaw::CheckRefresh() {
  if (D.UI[ForceLawPreset__].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawTailTol_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawNormali_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawScale___].hasChanged()) isRefreshed= false;
  for (int idxParam= ForceLawA_______; idxParam <= ForceLaw30______; idxParam++)
    if (D.UI[idxParam].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void ParticForceLaw::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Reset data arrays
  Ref.clear();
  Pos.clear();
  Vel.clear();
  Acc.clear();
  For.clear();
  Col.clear();
  Sensor.clear();
  BCPos.clear();
  BCVel.clear();
  BCFor.clear();
  SpatialSort.clear();

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
  if (D.Plot.size() < 4) D.Plot.resize(4);
  D.Plot[2].name= "KE";
  D.Plot[3].name= "ForceMag";
  if (D.Plot[2].val.size() < 10000) {
    D.Plot[2].val.reserve(10000);
    D.Plot[3].val.reserve(10000);
    D.Plot[2].val.push_back(0.0f);
    D.Plot[3].val.push_back(0.0f);
    for (int k= 0; k < (int)Pos.size(); k++) {
      D.Plot[2].val[D.Plot[2].val.size() - 1]+= ParticleMass * Vel[k].normSquared();
      D.Plot[3].val[D.Plot[3].val.size() - 1]+= For[k].norm();
    }
  }
}


// Draw the project
void ParticForceLaw::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Display particles
  if (D.displayMode1) {
    if (D.UI[VisuSimple______].B()) {
      glPointSize(1000.0f * D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F());
      glBegin(GL_POINTS);
      for (int k= 0; k < (int)Pos.size(); k++) {
        glColor3fv(Col[k].array());
        glVertex3fv(Pos[k].array());
      }
      glEnd();
    }
    else {
      glEnable(GL_LIGHTING);
      for (int k= 0; k < (int)Pos.size(); k++) {
        glColor3fv(Col[k].array());
        glPushMatrix();
        glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
        glutSolidSphere(D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F(), 12, 6);
        glPopMatrix();
      }
      glDisable(GL_LIGHTING);
    }
  }

  // Draw the spatial sort status
  if (D.displayMode2) {
    int nX, nY, nZ, nB;
    Field::GetFieldDimensions(SpatialSort, nX, nY, nZ, nB);
    glPointSize(1000.0f * D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F());
    glBegin(GL_POINTS);
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          float r= 0.5, g= 0.5, b= 0.5;
          Colormap::RatioToJetBrightSmooth((float)SpatialSort[x][y][z].size() / D.UI[SpatialSortSize_].F(), r, g, b);
          glColor3f(r, g, b);
          glVertex3f(D.boxMin[0] + (D.boxMax[0] - D.boxMin[0]) * ((float)x + 0.5f) / (float)nX,
                     D.boxMin[1] + (D.boxMax[1] - D.boxMin[1]) * ((float)y + 0.5f) / (float)nY,
                     D.boxMin[2] + (D.boxMax[2] - D.boxMin[2]) * ((float)z + 0.5f) / (float)nZ);
        }
      }
    }
    glEnd();
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("DrawT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Generate the base point cloud with an FCC pattern
  std::vector<Vec::Vec3<float>> points;
  const float minDist= std::sqrt(2.0f) / 2.0f;
  const int nX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * minDist));
  const int nY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * minDist));
  const int nZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * minDist));
  for (int x= 0; x < (nX / 2) * 2 + 1; x++) {
    for (int y= 0; y < (nY / 2) * 2 + 1; y++) {
      for (int z= 0; z < (nZ / 2) * 2 + 1; z++) {
        if (((x + y + z) % 2 == 0))
          points.push_back(Vec::Vec3<float>(
              D.boxMin[0] + x * D.UI[LatticePitch____].F() * minDist,
              D.boxMin[1] + y * D.UI[LatticePitch____].F() * minDist,
              D.boxMin[2] + z * D.UI[LatticePitch____].F() * minDist));
      }
    }
  }

  // Recenter the point cloud
  Vec::Vec3<float> avgPos(0.0f, 0.0f, 0.0f);
  for (int k= 0; k < (int)points.size(); k++)
    avgPos= avgPos + points[k];
  avgPos/= (float)points.size();
  for (int k= 0; k < (int)points.size(); k++)
    for (int dim= 0; dim < 3; dim++)
      points[k][dim]= points[k][dim] - avgPos[dim] + (D.boxMin[dim] + D.boxMax[dim]) / 2.0f;

  // Get the points in the box
  oPointCloud.clear();
  for (int k= 0; k < (int)points.size(); k++)
    if (points[k][0] >= D.boxMin[0] && points[k][0] <= D.boxMax[0])
      if (points[k][1] >= D.boxMin[1] && points[k][1] <= D.boxMax[1])
        if (points[k][2] >= D.boxMin[2] && points[k][2] <= D.boxMax[2])
          oPointCloud.push_back(points[k]);

  // Get the particle mass from the cloud and target density
  const float boxVolume= (D.boxMax[0] - D.boxMin[0]) * (D.boxMax[1] - D.boxMin[1]) * (D.boxMax[2] - D.boxMin[2]);
  ParticleMass= boxVolume * D.UI[MaterialDensity_].F() / (float)oPointCloud.size();

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
    if (D.UI[Scenario2DID____].I() == 1) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/CouponsSet.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 2) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/SandiaFracture.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 3) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/KalthoffFracture.bmp", imageRGBA, false);
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
            if (colRGBA[0] < 0.1f) Sensor.push_back(1);
            if (colRGBA[0] > 0.9f) BCPos.push_back(1);
            if (colRGBA[1] < 0.1f) BCVel.push_back(-1);
            if (colRGBA[1] > 0.9f) BCVel.push_back(1);
            if (colRGBA[2] < 0.1f) BCFor.push_back(-1);
            if (colRGBA[2] > 0.9f) BCFor.push_back(1);
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
        BCPos.push_back(1);
      }
      else if (RelPos[0] < 0.6f && RelPos[2] < 0.06f && RelPos[1] > 0.5f) {
        Pos.push_back(iPointCloud[k]);
        BCPos.push_back(1);
      }
    }
    // Ball blasting through wall
    else if (D.UI[ScenarioPreset__].I() == 4) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.2f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
      }
      else if (RelPos[2] > 0.60f && RelPos[2] < 0.65f) {
        Pos.push_back(iPointCloud[k]);
      }
    }
    // Balls colliding
    else if (D.UI[ScenarioPreset__].I() == 5) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.45f, 0.2f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
      }
      else if ((RelPos - Vec::Vec3<float>(0.5f, 0.55f, 0.8f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F()));
      }
    }
    // Coupon stretch - Velocity or Force driven
    else if (D.UI[ScenarioPreset__].I() == 6 || D.UI[ScenarioPreset__].I() == 7) {
      if (RelPos[0] > 0.45f && RelPos[0] < 0.55f) {
        if (RelPos[1] > 0.35f && RelPos[1] < 0.65f) {
          if (RelPos[2] > 0.2f && RelPos[2] < 0.8f) {
            Pos.push_back(iPointCloud[k]);
            if (D.UI[ScenarioPreset__].I() == 6) {
              if (RelPos[2] < 0.3f) {
                Sensor.push_back(1);
                BCVel.push_back(-1);
              }
              else if (RelPos[2] > 0.7f) {
                BCVel.push_back(1);
              }
            }
            else {
              if (RelPos[2] < 0.3f) {
                Sensor.push_back(1);
                BCFor.push_back(-1);
              }
              else if (RelPos[2] > 0.7f) {
                BCFor.push_back(1);
              }
            }
          }
        }
      }
    }
    // 3 Point flexture - Velocity or Force driven
    else if (D.UI[ScenarioPreset__].I() == 8 || D.UI[ScenarioPreset__].I() == 9) {
      if (RelPos[1] > 0.1f && RelPos[1] < 0.9f) {
        if ((RelPos - Vec::Vec3<float>(RelPos[0], 0.15f, 0.38f)).norm() < 0.05f ||
            (RelPos - Vec::Vec3<float>(RelPos[0], 0.85f, 0.38f)).norm() < 0.05f) {
          Pos.push_back(iPointCloud[k]);
          BCPos.push_back(1);
        }
        else if ((RelPos - Vec::Vec3<float>(RelPos[0], 0.5f, 0.62f)).norm() < 0.05f) {
          Pos.push_back(iPointCloud[k]);
          Sensor.push_back(1);
          if (D.UI[ScenarioPreset__].I() == 8) BCVel.push_back(-1);
          if (D.UI[ScenarioPreset__].I() == 9) BCFor.push_back(-1);
        }
        else if (RelPos[2] > 0.45f && RelPos[2] < 0.55f) {
          if (std::abs(RelPos[0] - 0.5f) < 0.5f * D.UI[Scenario2DThick_].F()) {
            Pos.push_back(iPointCloud[k]);
          }
        }
      }
    }

    // Fill the missing default values
    if (Sensor.size() < Pos.size()) Sensor.push_back(0);
    if (BCPos.size() < Pos.size()) BCPos.push_back(0);
    if (BCVel.size() < Pos.size()) BCVel.push_back(0);
    if (BCFor.size() < Pos.size()) BCFor.push_back(0);

    if (Ref.size() < Pos.size()) Ref.push_back(Pos[Pos.size() - 1]);
    if (Vel.size() < Pos.size()) Vel.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Acc.size() < Pos.size()) Acc.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (For.size() < Pos.size()) For.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Col.size() < Pos.size()) Col.push_back(Vec::Vec3<float>(0.5f, 0.5f, 0.5f));
  }

  SimTime= 0.0f;

  if (D.UI[VerboseLevel____].I() >= 1) printf("ScenarioT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildForceLaw() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the force law
  ForceLaw.clear();

  // Custom force law
  if (D.UI[ForceLawPreset__].I() == 0) {
    // Create the force law as two a smooth polylines
    std::vector<Vec::Vec3<double>> PolylineA;
    PolylineA.push_back(Vec::Vec3<double>{0.0, 0.0, D.UI[ForceLawA_______].D()});
    PolylineA.push_back(Vec::Vec3<double>{0.0, 0.8, D.UI[ForceLaw08______].D()});
    PolylineA.push_back(Vec::Vec3<double>{0.0, 0.9, D.UI[ForceLaw09______].D()});
    PolylineA.push_back(Vec::Vec3<double>{0.0, 1.0, D.UI[ForceLawB_______].D()});
    std::vector<Vec::Vec3<double>> PolylineB;
    PolylineB.push_back(Vec::Vec3<double>{0.0, 1.0, D.UI[ForceLawB_______].D()});
    PolylineB.push_back(Vec::Vec3<double>{0.0, 1.1, D.UI[ForceLaw11______].D()});
    PolylineB.push_back(Vec::Vec3<double>{0.0, 1.2, D.UI[ForceLaw12______].D()});
    PolylineB.push_back(Vec::Vec3<double>{0.0, 1.3, D.UI[ForceLaw13______].D()});
    PolylineB.push_back(Vec::Vec3<double>{0.0, std::sqrt(2), D.UI[ForceLawC_______].D()});
    std::vector<Vec::Vec3<double>> PolylineC;
    PolylineC.push_back(Vec::Vec3<double>{0.0, std::sqrt(2), D.UI[ForceLawC_______].D()});
    PolylineC.push_back(Vec::Vec3<double>{0.0, 1.5, D.UI[ForceLaw15______].D()});
    PolylineC.push_back(Vec::Vec3<double>{0.0, 2.0, D.UI[ForceLaw20______].D()});
    PolylineC.push_back(Vec::Vec3<double>{0.0, 2.5, D.UI[ForceLaw25______].D()});
    PolylineC.push_back(Vec::Vec3<double>{0.0, 3.0, D.UI[ForceLaw30______].D()});
    Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineA);
    Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineB);
    Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineC);
    std::vector<Vec::Vec3<double>> Polyline;
    Polyline.insert(Polyline.end(), PolylineA.begin(), PolylineA.end());
    Polyline.insert(Polyline.end(), PolylineB.begin(), PolylineB.end());
    Polyline.insert(Polyline.end(), PolylineC.begin(), PolylineC.end());

    // Sample the force value at fixed distance intervals
    ForceLawStep= std::sqrt(2.0f) / 128.0f;
    ForceLaw.resize((int)std::ceil(3.0f / ForceLawStep));
    for (int k= 0; k < (int)ForceLaw.size(); k++) {
      float dist= ForceLawStep * float(k);
      int idxLow= 0, idxUpp= 0;
      for (int idxVert= 0; idxVert < (int)Polyline.size() - 1; idxVert++) {
        idxLow= idxVert;
        idxUpp= std::min(idxVert + 1, (int)Polyline.size() - 1);
        if (Polyline[idxLow][1] <= dist && Polyline[idxUpp][1] >= dist)
          break;
      }
      double ratio= 0.0;
      if (Polyline[idxUpp][1] - Polyline[idxLow][1] > 0.0)
        ratio= (dist - Polyline[idxLow][1]) / (Polyline[idxUpp][1] - Polyline[idxLow][1]);
      ForceLaw[k]= (1.0 - ratio) * Polyline[idxLow][2] + ratio * Polyline[idxUpp][2];
    }
  }
  else {
    ForceLawStep= 0.05f;
    //                             00      05      10      15      20      25      30      35      40      45      50      55      60      65      70      75      80      85      90      95      0       05       10       15       20       25       30       35       40       45       50       55       60       65       70       75       80       85       90       95       00       05       10       15       20       25       30       35       40       45       50       55       60       65       70       75       80       85       90       95      00
    if (D.UI[ForceLawPreset__].I() == 2)  // Elastic (fig 2 MIT paper)
      ForceLaw= std::vector<float>{1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 0.5e6f, 0.0e6f, -0.5e6f, -1.0e6f, -0.5e6f, 0.0e6f, 0.5e6f, 1.0e6f, 0.5e6f, 0.0e6f};
    else if (D.UI[ForceLawPreset__].I() == 3)  // Brittle (fig 3 MIT paper)
      ForceLaw= std::vector<float>{1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 1.0e6f, 0.5e6f, 0.0e6f, -0.1e6f, 0.0e6f, 0.07e6f, 0.10e6f, 0.12e6f, 0.08e6f, 0.03e6f, 0.0e6f};
    else if (D.UI[ForceLawPreset__].I() == 4)  // Viscous (fig 4 MIT paper)
      ForceLaw= std::vector<float>{100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 95.0, 88.0, 70.0, 50.0, 25.0, 0.0, -20.0, -28.0, -30.0, -28.0, -25.0, -20.0, -12.0, -7.0, -5.0, 0.0};
    else if (D.UI[ForceLawPreset__].I() == 5)  // Steel AISI 4340 (fig 9 MIT paper)
      ForceLaw= std::vector<float>{3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 3.0e10f, 2.0e10f, 0.0e10f, -2.0e10f, -2.5e10f, -2.0e10f, 0.0e10f, 2.0e10f, 2.5e10f, 2.0e10f, 0.0e10f};
    else if (D.UI[ForceLawPreset__].I() == 6)  // Delrin
      ForceLaw= std::vector<float>{88.e3f, 88.e3f, 88.e3f, 87.e3f, 87.e3f, 86.e3f, 85.e3f, 85.e3f, 83.e3f, 81.e3f, 78.e3f, 75.e3f, 72.e3f, 66.e3f, 60.e3f, 53.e3f, 45.e3f, 37.e3f, 29.e3f, 20.e3f, 0.e3f, -30.e3f, -32.e3f, -33.e3f, -33.e3f, -33.e3f, -32.e3f, -32.e3f, -31.e3f, -31.e3f, -30.e3f, -30.e3f, -30.e3f, -29.e3f, -29.e3f, -28.e3f, -28.e3f, -28.e3f, -28.e3f, -29.e3f, -29.e3f, -30.e3f, -30.e3f, -30.e3f, -31.e3f, -32.e3f, -33.e3f, -35.e3f, -36.e3f, -38.e3f, -40.e3f, -42.e3f, -43.e3f, -44.e3f, -44.e3f, -43.e3f, -41.e3f, -38.e3f, -32.e3f, -21.e3f, 0.e3f};
    else  // Default force law (fig 1 MIT paper)
      ForceLaw= std::vector<float>{1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.5e8f, 1.0e8f, 0.0e8f, -1.5e8f, -2.3e8f, -1.5e8f, 0.0e8f, 1.1e8f, 1.0e8f, 0.5e8f, 0.0e8f};
  }

  // Cutoff the zero tail of the force law
  if (D.UI[ForceLawTailTol_].F() > 0.0f) {
    float BaseVal= ForceLaw[0];
    for (int k= (int)ForceLaw.size() - 2; k >= 0; k--) {
      if (std::abs(ForceLaw[k + 1]) < std::abs(D.UI[ForceLawTailTol_].F() * BaseVal) &&
          std::abs(ForceLaw[k]) < std::abs(D.UI[ForceLawTailTol_].F() * BaseVal))
        ForceLaw.pop_back();
      else
        break;
    }
  }

  // Optionally normalize the force law
  if (D.UI[ForceLawNormali_].B()) {
    float BaseVal= ForceLaw[0];
    for (int k= 0; k < (int)ForceLaw.size(); k++)
      ForceLaw[k]/= BaseVal;
  }

  // Scale the force law
  for (int k= 0; k < (int)ForceLaw.size(); k++)
    ForceLaw[k]*= D.UI[ForceLawScale___].F();

  if (D.UI[VerboseLevel____].I() >= 1) printf("ForceLawT %f\n", Timer::PopTimer());
}


void ParticForceLaw::ComputeSpatialSort() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Check validity of inputs
  if (D.UI[SpatialSortNb___].I() < 1) return;
  if (D.UI[SpatialSortSize_].I() < 1) return;

  // Compute the appropriate voxel size and grid resolution
  const float boxVolume= (D.boxMax[0] - D.boxMin[0]) * (D.boxMax[1] - D.boxMin[1]) * (D.boxMax[2] - D.boxMin[2]);
  const float stepSize= std::pow(boxVolume, 1.0f / 3.0f) / std::pow((float)D.UI[SpatialSortNb___].I(), 1.0f / 3.0f);
  const int nX= std::max(1, (int)std::round((D.boxMax[0] - D.boxMin[0]) / stepSize));
  const int nY= std::max(1, (int)std::round((D.boxMax[1] - D.boxMin[1]) / stepSize));
  const int nZ= std::max(1, (int)std::round((D.boxMax[2] - D.boxMin[2]) / stepSize));

  // Allocate/initialize the spatial sort if needed
  int nXOld, nYOld, nZOld, nBOld;
  Field::GetFieldDimensions(SpatialSort, nXOld, nYOld, nZOld, nBOld);
  if (nXOld != nX || nYOld != nY || nZOld != nZ || nBOld != D.UI[SpatialSortSize_].I())
    SpatialSort= Field::AllocField4D(nX, nY, nZ, D.UI[SpatialSortSize_].I(), -1);

  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        SpatialSort[x][y][z].clear();

  // Fill the spatial sort
  for (int k= 0; k < (int)Pos.size(); k++) {
    const int idxX= (int)std::floor((float)nX * (Pos[k][0] - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxY= (int)std::floor((float)nY * (Pos[k][1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZ= (int)std::floor((float)nZ * (Pos[k][2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    if (idxX >= 0 && idxX < nX && idxY >= 0 && idxY < nY && idxZ >= 0 && idxZ < nZ)
      if ((int)SpatialSort[idxX][idxY][idxZ].size() < D.UI[SpatialSortSize_].I())
        SpatialSort[idxX][idxY][idxZ].push_back(k);
  }

  // Plot spatial sort occupancy
  if (D.Plot.size() < 1) D.Plot.resize(1);
  D.Plot[0].name= "SpatialSort";
  D.Plot[0].val.resize(nX * nY * nZ + 2);
  D.Plot[0].val[nX * nY * nZ + 0]= 0;
  D.Plot[0].val[nX * nY * nZ + 1]= D.UI[SpatialSortSize_].I();
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        D.Plot[0].val[x * nY * nZ + y * nZ + z]= (double)SpatialSort[x][y][z].size();

  if (D.UI[VerboseLevel____].I() >= 1) printf("SpatialSortT %f\n", Timer::PopTimer());
}


void ParticForceLaw::ComputeForces() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Check validity of inputs
  if ((int)ForceLaw.size() <= 1) return;

  // Reset forces
  std::fill(For.begin(), For.end(), Vec::Vec3<float>(0.0, 0.0, 0.0));

  // Precompute values
  const float forceReach= D.UI[LatticePitch____].F() * ForceLawStep * (float)(ForceLaw.size() - 1);

  // Compute the spatial sort for linear neighbor search
  ComputeSpatialSort();
  int nX, nY, nZ, nB;
  Field::GetFieldDimensions(SpatialSort, nX, nY, nZ, nB);

// Compute particle forces
#pragma omp parallel for
  for (int k0= 0; k0 < (int)Pos.size(); k0++) {
    // Get range to check in spatial sort
    const int idxXBeg= (int)std::floor((float)nX * (Pos[k0][0] - forceReach - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxYBeg= (int)std::floor((float)nY * (Pos[k0][1] - forceReach - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZBeg= (int)std::floor((float)nZ * (Pos[k0][2] - forceReach - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    const int idxXEnd= (int)std::floor((float)nX * (Pos[k0][0] + forceReach - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxYEnd= (int)std::floor((float)nY * (Pos[k0][1] + forceReach - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZEnd= (int)std::floor((float)nZ * (Pos[k0][2] + forceReach - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    // Check range in spatial sort
    for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, nX - 1); x++) {
      for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, nY - 1); y++) {
        for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, nZ - 1); z++) {
          // Check candidate particles
          for (int k1 : SpatialSort[x][y][z]) {
            if (k1 == k0) continue;
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
            // Apply inter-particle damping linearly proportional to difference in radial velocity of particle pair
            // TODO make damping distance dependant ?
            const float RadialVel= (Vel[k0] - Vel[k1]).dot(distVec / distVal);
            For[k0]-= D.UI[VelocityDamping_].F() * RadialVel * distVec / distVal;
          }
        }
      }
    }
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("ForcesT %f\n", Timer::PopTimer());
}


void ParticForceLaw::ApplyBCForces() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Precompute values
  const float dt= D.UI[TimeStep________].F();
  const Vec::Vec3<float> BCForVecNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCForVecPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCVelVecNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> BCVelVecPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());

  // Apply boundary conditions via force controller
  for (int k= 0; k < (int)Pos.size(); k++) {
    if (BCPos[k] != 0) {
      // xt+1 = xt + dt*vt + dt*dt*ft/m
      // xt+1 == xbc  =>  xt + dt*vt + dt*dt*ft/m == xbc
      //              =>  dt*dt*ft/m              == xbc - xt - dt*vt
      //              =>  ft                      == m*(xbc - xt - dt*vt)/dt*dt
      Vec::Vec3<float> ErrVec= Ref[k] - Pos[k];
      For[k]+= D.UI[BCPosCoeff______].F() * ParticleMass * (ErrVec.normalized() * ErrVec.normSquared() - dt * Vel[k]) / (dt * dt);
    }
    else if (BCVel[k] != 0) {
      // vt+1 = vt + dt*ft/m
      // vt+1 == vbc  =>  vt + dt*ft/m == vbc
      //              =>  dt*ft/m      == vbc - vt
      //              =>  ft           == m*(vbc-vt)/dt
      Vec::Vec3<float> ErrVec= ((BCVel[k] < 0) ? (BCVelVecNega) : (BCVelVecPosi)) - Vel[k];
      For[k]+= D.UI[BCVelCoeff______].F() * ParticleMass * ErrVec.normalized() * ErrVec.normSquared() / dt;
    }
    else if (BCFor[k] != 0) {
      For[k]+= (BCFor[k] < 0) ? (BCForVecNega) : (BCForVecPosi);
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("BoundCondT %f\n", Timer::PopTimer());
}


void ParticForceLaw::StepSimulation() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Precompute values
  const float dt= D.UI[TimeStep________].F();
  const Vec::Vec3<float> BCForVecNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCForVecPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCVelVecNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> BCVelVecPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());

  // Explicit integration
  if (D.UI[IntegType_______].I() == 0) {
    // Evaluate net forces acting on particles
    ComputeForces();                                                                   // ft
    if (D.UI[UseForceControl_].B()) ApplyBCForces();                                   // Boundary conditions with force controller
    else                                                                               // Boundary conditions check
      for (int k= 0; k < (int)Pos.size(); k++)                                         // Loop through elements
        if (BCFor[k] != 0) For[k]+= (BCFor[k] < 0) ? (BCForVecNega) : (BCForVecPosi);  // Add force
    // Euler integration
    for (int k= 0; k < (int)Pos.size(); k++) {
      Acc[k]= For[k] / ParticleMass;                           // at = ft / m
      Vel[k]+= dt * Acc[k];                                    // vt+1 = vt + dt * at
      if (!D.UI[UseForceControl_].B() && BCVel[k] != 0)        // Boundary conditions check
        Vel[k]= (BCVel[k] < 0) ? BCVelVecNega : BCVelVecPosi;  // Overwrite velocity
      Pos[k]+= dt * Vel[k];                                    // xt+1 = xt + dt * vt+1
      if (!D.UI[UseForceControl_].B() && BCPos[k] != 0) {      // Boundary conditions check
        Pos[k]= Ref[k];                                        // Overwrite Position
        Vel[k]= Vec::Vec3<float>{0.0f, 0.0f, 0.0f};            // Overwrite Velocity
      }
    }
  }
  else {
    // Velocity Verlet integration - position update
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (!D.UI[UseForceControl_].B() && BCPos[k] != 0) {  // Boundary conditions check
        Pos[k]= Ref[k];                                    // Overwrite Position
        Vel[k]= Vec::Vec3<float>{0.0f, 0.0f, 0.0f};        // Overwrite Velocity
      }
      else {
        Pos[k]+= dt * Vel[k] + 0.5f * dt * dt * Acc[k];  // xt+1 = xt + dt * vt + 0.5 * dt * dt * at
      }
    }
    // Evaluate net forces acting on particles
    ComputeForces();                                                                   // ft+1
    if (D.UI[UseForceControl_].B()) ApplyBCForces();                                   // Boundary conditions with force controller
    else                                                                               // Boundary conditions check
      for (int k= 0; k < (int)Pos.size(); k++)                                         // Loop through elements
        if (BCFor[k] != 0) For[k]+= (BCFor[k] < 0) ? (BCForVecNega) : (BCForVecPosi);  // Add force
    // Velocity Verlet integration - acceleration and velocity update
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (!D.UI[UseForceControl_].B() && BCVel[k] != 0) {      // Boundary conditions check
        Vel[k]= (BCVel[k] < 0) ? BCVelVecNega : BCVelVecPosi;  // Overwrite velocity
        Acc[k]= 0.0f;                                          // Overwrite acceleration
      }
      else {
        const Vec::Vec3<float> oldAcc= Acc[k];   // Store previous acceleration
        Acc[k]= For[k] / ParticleMass;           // at+1 = ft+1 / m
        Vel[k]+= 0.5f * dt * (oldAcc + Acc[k]);  // vt+1 = vt + 0.5 * dt * (at + at+1)
      }
    }
  }

  // Set particle color with gradual decay
  for (int k= 0; k < (int)Pos.size(); k++) {
    float r= 0.5f, g= 0.5f, b= 0.5f;
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
    if (D.UI[ColorMode_______].I() == 3) Colormap::RatioToJetBrightSmooth(For[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
    Col[k]= (1.0f - D.UI[ColorDecay______].F()) * Col[k] + D.UI[ColorDecay______].F() * Vec::Vec3<float>(r, g, b);
  }

  // Write the status
  SimTime+= dt;
  if ((int)D.Status.size() < 3) D.Status.resize(3);
  D.Status[0]= std::string{"NbParticles: "} + std::to_string((int)Pos.size());
  D.Status[1]= std::string{"ParticleMass: "} + std::to_string(ParticleMass);
  D.Status[2]= std::string{"SimTime: "} + std::to_string(SimTime);

  // Scatter plot of sensor data
  if (D.Scatter.size() < 1) D.Scatter.resize(1);
  D.Scatter[0].name= "ForceDisp";
  float sumPos= 0.0f;
  float sumFor= 0.0f;
  int count= 0;
  for (int k= 0; k < (int)Pos.size(); k++) {
    if (Sensor[k]) {
      count++;
      sumPos+= Pos[k][2];
      sumFor+= For[k].norm();
    }
  }
  D.Scatter[0].val.push_back(std::array<double, 2>{sumPos / (float)count, sumFor / (float)count});

  if (D.UI[VerboseLevel____].I() >= 1) printf("StepT %f\n", Timer::PopTimer());
}
