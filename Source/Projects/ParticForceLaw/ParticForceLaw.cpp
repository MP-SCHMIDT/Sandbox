#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <numbers>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/BoxGrid.hpp"
#include "Geom/MarchingCubes.hpp"
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
    D.UI.push_back(ParamUI("LatticePattern__", 2));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("ForceLawPreset0_", 0));
    D.UI.push_back(ParamUI("ForceLawPreset1_", 0));
    D.UI.push_back(ParamUI("ForceLawPreset2_", 0));
    D.UI.push_back(ParamUI("ForceLawScale___", 100.0));
    D.UI.push_back(ParamUI("ForceLawA_______", 1.0));
    D.UI.push_back(ParamUI("ForceLaw08______", 1.0));
    D.UI.push_back(ParamUI("ForceLaw09______", 1.0));
    D.UI.push_back(ParamUI("ForceLawB_______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw11______", -1.0));
    D.UI.push_back(ParamUI("ForceLaw12______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw13______", 1.0));
    D.UI.push_back(ParamUI("ForceLawC_______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw15______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw20______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw25______", 0.0));
    D.UI.push_back(ParamUI("ForceLaw30______", 0.0));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("BCVelX__________", 0.0));
    D.UI.push_back(ParamUI("BCVelY__________", 0.0));
    D.UI.push_back(ParamUI("BCVelZ__________", 1.0));
    D.UI.push_back(ParamUI("BCForX__________", 0.0));
    D.UI.push_back(ParamUI("BCForY__________", 0.0));
    D.UI.push_back(ParamUI("BCForZ__________", 1.0));
    D.UI.push_back(ParamUI("StepsPerDraw____", 1));
    D.UI.push_back(ParamUI("TimeStep________", 0.002));
    D.UI.push_back(ParamUI("MaterialDensity_", 200.0));
    D.UI.push_back(ParamUI("DampingRadRel___", 0.1));
    D.UI.push_back(ParamUI("DampingVelRel___", 0.0));
    D.UI.push_back(ParamUI("BucketsCount____", 1000));
    D.UI.push_back(ParamUI("BucketsCapacity_", 40));
    D.UI.push_back(ParamUI("IntegType_______", 1));
    D.UI.push_back(ParamUI("UseForceControl_", 0));
    D.UI.push_back(ParamUI("BCPosCoeff______", 1.0));
    D.UI.push_back(ParamUI("BCVelCoeff______", 1.0));
    D.UI.push_back(ParamUI("MetaballVoxSize_", 0.04));
    D.UI.push_back(ParamUI("MetaballIsoval__", 0.5));
    D.UI.push_back(ParamUI("ColorMode_______", 3));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("ColorDecay______", 0.5));
    D.UI.push_back(ParamUI("VisuScale_______", 0.5));
    D.UI.push_back(ParamUI("VisuSimple______", 1));
    D.UI.push_back(ParamUI("VisuShowOOB_____", 0));
    D.UI.push_back(ParamUI("TestParam0______", 0));
    D.UI.push_back(ParamUI("TestParam1______", 0));
    D.UI.push_back(ParamUI("TestParam2______", 0));
    D.UI.push_back(ParamUI("TestParam3______", 0));
    D.UI.push_back(ParamUI("TestParam4______", 0));
    D.UI.push_back(ParamUI("TestParam5______", 0));
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
  if (D.UI[LatticePattern__].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool ParticForceLaw::CheckRefresh() {
  if (D.UI[ForceLawPreset0_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawPreset1_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawPreset2_].hasChanged()) isRefreshed= false;
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
  Mat.clear();
  Sensor.clear();
  BCPos.clear();
  BCVel.clear();
  BCFor.clear();
  Buckets.clear();
  nX= nY= nZ= 0;
  MetaballIsUpdated= false;
  Verts.clear();
  Tris.clear();

  // Get domain dimensions
  D.boxMin= {0.5 - 0.5 * D.UI[DomainX_________].D(), 0.5 - 0.5 * D.UI[DomainY_________].D(), 0.5 - 0.5 * D.UI[DomainZ_________].D()};
  D.boxMax= {0.5 + 0.5 * D.UI[DomainX_________].D(), 0.5 + 0.5 * D.UI[DomainY_________].D(), 0.5 + 0.5 * D.UI[DomainZ_________].D()};

  // Generate the full point cloud over the domain
  std::vector<Vec::Vec3<float>> pointCloud;
  BuildBaseCloud(pointCloud);

  // Generate the scenario
  BuildScenario(pointCloud);

  // Generate the spatial partition
  ComputeBuckets();
}


// Refresh the project
void ParticForceLaw::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Generate the force law
  BuildForceLaws();

  // Draw the force law in the plot
  if (D.Plot.size() < 6) D.Plot.resize(6);
  for (int idxMat= 0; idxMat < (int)ForceLaws.size(); idxMat++) {
    D.Plot[3 + idxMat].name= "ForceLaw";
    D.Plot[3 + idxMat].val.resize(ForceLaws[idxMat].size());
    for (int k= 0; k < (int)ForceLaws[idxMat].size(); k++)
      D.Plot[3 + idxMat].val[k]= ForceLaws[idxMat][k];
  }
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
    MetaballIsUpdated= false;
  }

  // Plot data
  if (D.Plot.size() < 3) D.Plot.resize(3);
  D.Plot[1].name= "KE";
  D.Plot[2].name= "ForceMag";
  if (D.Plot[1].val.size() < 10000) {
    D.Plot[1].val.reserve(10000);
    D.Plot[2].val.reserve(10000);
    D.Plot[1].val.push_back(0.0f);
    D.Plot[2].val.push_back(0.0f);
    for (int k= 0; k < (int)Pos.size(); k++) {
      D.Plot[1].val[D.Plot[1].val.size() - 1]+= Vel[k].normSquared();
      D.Plot[2].val[D.Plot[2].val.size() - 1]+= For[k].norm();
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
        if (!D.UI[VisuShowOOB_____].B())
          if (Pos[k][0] < D.boxMin[0] || Pos[k][0] > D.boxMax[0] ||
              Pos[k][1] < D.boxMin[1] || Pos[k][1] > D.boxMax[1] ||
              Pos[k][2] < D.boxMin[2] || Pos[k][2] > D.boxMax[2]) continue;
        glColor3fv(Col[k].array());
        glVertex3fv(Pos[k].array());
      }
      glEnd();
    }
    else {
      glEnable(GL_LIGHTING);
      for (int k= 0; k < (int)Pos.size(); k++) {
        if (!D.UI[VisuShowOOB_____].B())
          if (Pos[k][0] < D.boxMin[0] || Pos[k][0] > D.boxMax[0] ||
              Pos[k][1] < D.boxMin[1] || Pos[k][1] > D.boxMax[1] ||
              Pos[k][2] < D.boxMin[2] || Pos[k][2] > D.boxMax[2]) continue;
        glColor3fv(Col[k].array());
        glPushMatrix();
        glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
        glutSolidSphere(D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F(), 12, 6);
        glPopMatrix();
      }
      glDisable(GL_LIGHTING);
    }
  }

  // Display spatial partition buckets status
  if (!D.displayMode2) {
    glLineWidth(2.0f);
    // Get dimensions
    double stepX, stepY, stepZ;
    BoxGrid::GetVoxelSizes(nX, nY, nZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ);
    const int bucketCapacity= std::max(D.UI[BucketsCapacity_].I(), 1);
    // Set transformation
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5f * (float)stepX, D.boxMin[1] + 0.5f * (float)stepY, D.boxMin[2] + 0.5f * (float)stepZ);
    glScalef((float)stepX, (float)stepY, (float)stepZ);
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          // Color by occupancy
          float r= 0.5f, g= 0.5f, b= 0.5f;
          Colormap::RatioToJetBrightSmooth((float)Buckets[x][y][z].size() / (float)bucketCapacity, r, g, b);
          glColor3f(r, g, b);
          // Draw wire box
          glPushMatrix();
          glTranslatef((float)x, (float)y, (float)z);
          glutWireCube(1.0);
          glPopMatrix();
        }
      }
    }
    glPopMatrix();
    glLineWidth(1.0f);
  }


  // Draw triangles
  if (!D.displayMode3) {
    if (!MetaballIsUpdated)
      ComputeMetaballs();
    glEnable(GL_LIGHTING);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_TRIANGLES);
    for (int k= 0; k < (int)Tris.size(); k++) {
      Vec::Vec3 v0(Verts[Tris[k][0]][0], Verts[Tris[k][0]][1], Verts[Tris[k][0]][2]);
      Vec::Vec3 v1(Verts[Tris[k][1]][0], Verts[Tris[k][1]][1], Verts[Tris[k][1]][2]);
      Vec::Vec3 v2(Verts[Tris[k][2]][0], Verts[Tris[k][2]][1], Verts[Tris[k][2]][2]);
      Vec::Vec3 n0= (v1 - v0).cross(v2 - v0).normalized();
      Vec::Vec3 n1= n0;
      Vec::Vec3 n2= n0;
      glNormal3f(n0[0], n0[1], n0[2]);
      glVertex3f(v0[0], v0[1], v0[2]);
      glNormal3f(n1[0], n1[1], n1[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
      glNormal3f(n2[0], n2[1], n2[2]);
      glVertex3f(v2[0], v2[1], v2[2]);
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }
  else {
    MetaballIsUpdated= false;
    Verts.clear();
    Tris.clear();
  }

  // Write the status
  if ((int)D.Status.size() < 3) D.Status.resize(3);
  D.Status[0]= std::string{" NbBuckets="} + std::to_string(nX * nY * nZ);
  D.Status[1]= std::string{" NbParticles="} + std::to_string((int)Pos.size());
  D.Status[2]= std::string{" SimTime="} + std::to_string(SimTime);

  if (D.UI[VerboseLevel____].I() >= 1) printf("DrawT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the base point cloud
  oPointCloud.clear();

  // Check parameters
  if (D.UI[LatticePitch____].F() <= 0.0f) return;

  // Regular lattice pattern
  if (D.UI[LatticePattern__].I() == 0 || D.UI[LatticePattern__].I() == 1 || D.UI[LatticePattern__].I() == 2) {
    float minDist= 0.0f;
    if (D.UI[LatticePattern__].I() == 0) minDist= 1.0f;                                    // SCC pattern
    if (D.UI[LatticePattern__].I() == 1) minDist= (2.0f / 3.0f) * std::sqrt(3.0f) / 2.0f;  // BCC pattern
    if (D.UI[LatticePattern__].I() == 2) minDist= std::sqrt(2.0f) / 2.0f;                  // FCC pattern
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * minDist));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * minDist));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * minDist));
    for (int x= 0; x < (cellNbX / 2) * 2 + 1; x++) {
      for (int y= 0; y < (cellNbY / 2) * 2 + 1; y++) {
        for (int z= 0; z < (cellNbZ / 2) * 2 + 1; z++) {
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
  }
  // Poisson sphere sampling
  else if (D.UI[LatticePattern__].I() == 3) {
    ComputeBuckets();
    const int bucketCapacity= std::max(D.UI[BucketsCapacity_].I(), 1);
    int failStreak= 0;
    const int maxAttempts= 1000;
    const float relMinDist= 0.9;
    while (failStreak < maxAttempts) {
      const Vec::Vec3<float> candidate(Random::Val(D.boxMin[0], D.boxMax[0]), Random::Val(D.boxMin[1], D.boxMax[1]), Random::Val(D.boxMin[2], D.boxMax[2]));
      bool keep= true;
      // Get range to check in spatial partition
      int idxXBeg, idxYBeg, idxZBeg, idxXEnd, idxYEnd, idxZEnd;
      Vec::Vec3<float> vecOffset(D.UI[LatticePitch____].F(), D.UI[LatticePitch____].F(), D.UI[LatticePitch____].F());
      GetBucketIdx(candidate - vecOffset, idxXBeg, idxYBeg, idxZBeg);
      GetBucketIdx(candidate + vecOffset, idxXEnd, idxYEnd, idxZEnd);
      // Check range in spatial partition
      for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, nX - 1) && keep; x++) {
        for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, nY - 1) && keep; y++) {
          for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, nZ - 1) && keep; z++) {
            // Check candidate particles
            for (int k : Buckets[x][y][z]) {
              // Skip if invalid neighbor
              if ((candidate - oPointCloud[k]).normSquared() < std::pow(D.UI[LatticePitch____].F() * relMinDist, 2.0f)) {
                keep= false;
                break;
              }
            }
          }
        }
      }
      if (keep) {
        failStreak= 0;
        oPointCloud.push_back(candidate);
        int idxX, idxY, idxZ;
        GetBucketIdx(candidate, idxX, idxY, idxZ);
        if (idxX >= 0 && idxX < nX && idxY >= 0 && idxY < nY && idxZ >= 0 && idxZ < nZ)
          if ((int)Buckets[idxX][idxY][idxZ].size() < bucketCapacity)
            Buckets[idxX][idxY][idxZ].push_back((int)oPointCloud.size() - 1);
      }
      else {
        failStreak++;
      }
    }
  }

  // Recenter the point cloud
  Vec::Vec3<float> avgPos(0.0f, 0.0f, 0.0f);
  for (int k= 0; k < (int)oPointCloud.size(); k++)
    avgPos= avgPos + oPointCloud[k];
  avgPos/= (float)oPointCloud.size();
  for (int k= 0; k < (int)oPointCloud.size(); k++)
    for (int dim= 0; dim < 3; dim++)
      oPointCloud[k][dim]= oPointCloud[k][dim] - avgPos[dim] + (D.boxMin[dim] + D.boxMax[dim]) / 2.0f;

  if (D.UI[VerboseLevel____].I() >= 1) printf("BaseCloudT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildScenario(const std::vector<Vec::Vec3<float>>& iPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Check parameters
  if (D.boxMax[0] <= D.boxMin[0]) return;
  if (D.boxMax[1] <= D.boxMin[1]) return;
  if (D.boxMax[2] <= D.boxMin[2]) return;

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
    if (D.UI[Scenario2DID____].I() == 2) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/KalthoffFracture.bmp", imageRGBA, false);
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
        Mat.push_back(0);
      }
      else if (RelPos[0] > 0.1f && RelPos[0] < 0.9f &&
               RelPos[1] > 0.1f && RelPos[1] < 0.9f &&
               RelPos[2] > 0.45f && RelPos[2] < 0.55f) {
        Pos.push_back(iPointCloud[k]);
        Mat.push_back(1);
      }
    }
    // Balls colliding
    else if (D.UI[ScenarioPreset__].I() == 5) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.45f, 0.2f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Mat.push_back(0);
      }
      else if ((RelPos - Vec::Vec3<float>(0.5f, 0.55f, 0.8f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F()));
        Mat.push_back(1);
      }
    }
    // Coupon stretch - Velocity or Force driven
    else if (D.UI[ScenarioPreset__].I() == 6 || D.UI[ScenarioPreset__].I() == 7) {
      if (RelPos[0] > 0.1f && RelPos[0] < 0.9f) {
        if (RelPos[1] > 0.1f && RelPos[1] < 0.9f) {
          if (RelPos[2] > 0.2f && RelPos[2] < 0.8f) {
            Pos.push_back(iPointCloud[k]);
            if (D.UI[ScenarioPreset__].I() == 6) {
              if (RelPos[2] < 0.3f) {
                BCVel.push_back(-1);
              }
              else if (RelPos[2] > 0.7f) {
                Sensor.push_back(1);
                BCVel.push_back(1);
              }
            }
            else {
              if (RelPos[2] < 0.3f) {
                BCFor.push_back(-1);
              }
              else if (RelPos[2] > 0.7f) {
                Sensor.push_back(1);
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
    // Simple sphere
    else if (D.UI[ScenarioPreset__].I() == 10) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.5f)).norm() < 0.4f) {
        Pos.push_back(iPointCloud[k]);
        BCFor.push_back(-1);
        Sensor.push_back(1);
      }
    }

    // Fill the missing default values
    if (Mat.size() < Pos.size()) Mat.push_back(0);
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


void ParticForceLaw::BuildForceLaws() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the force law
  ForceLaws.clear();
  ForceLaws.resize(3);
  ForceLawSteps.clear();
  ForceLawSteps.resize(3);
  ForceLawRanges.clear();
  ForceLawRanges.resize(3);

  for (int idxMat= 0; idxMat < (int)ForceLaws.size(); idxMat++) {
    int presetID= -1;
    if (idxMat == 0) presetID= D.UI[ForceLawPreset0_].I();
    if (idxMat == 1) presetID= D.UI[ForceLawPreset1_].I();
    if (idxMat == 2) presetID= D.UI[ForceLawPreset2_].I();
    // Custom force law
    if (presetID == 0) {
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
      ForceLawSteps[idxMat]= std::sqrt(2.0f) / 128.0f;
      ForceLaws[idxMat].resize((int)std::ceil(3.0f / ForceLawSteps[idxMat]));
      for (int k= 0; k < (int)ForceLaws[idxMat].size(); k++) {
        float dist= ForceLawSteps[idxMat] * float(k);
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
        ForceLaws[idxMat][k]= (1.0 - ratio) * Polyline[idxLow][2] + ratio * Polyline[idxUpp][2];
      }
    }
    else {
      ForceLawSteps[idxMat]= 0.05f;
      ForceLaws[idxMat]= std::vector<float>{1.0f};
      // Hard coded force laws from MIT paper
      if (presetID == 1) ForceLaws[idxMat]= std::vector<float>{1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +5.00E+05, +0.00E+00, -5.00E+05, -1.00E+06, -5.00E+05, +0.00E+00, +5.00E+05, +1.00E+06, +5.00E+05, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};
      if (presetID == 2) ForceLaws[idxMat]= std::vector<float>{1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +5.00E+05, +0.00E+00, -1.00E+05, +0.00E+00, +7.00E+04, +1.00E+05, +1.20E+05, +8.00E+04, +3.00E+04, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};
      if (presetID == 3) ForceLaws[idxMat]= std::vector<float>{1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +9.50E+01, +8.80E+01, +7.00E+01, +5.00E+01, +2.50E+01, +0.00E+00, -2.00E+01, -2.80E+01, -3.00E+01, -2.80E+01, -2.50E+01, -2.00E+01, -1.20E+01, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};
      if (presetID == 4) ForceLaws[idxMat]= std::vector<float>{3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +2.00E+10, +0.00E+00, -2.00E+10, -2.50E+10, -2.00E+10, +0.00E+00, +2.00E+10, +2.50E+10, +2.00E+10, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};
      if (presetID == 5) ForceLaws[idxMat]= std::vector<float>{8.80E+04, +8.80E+04, +8.80E+04, +8.70E+04, +8.70E+04, +8.60E+04, +8.50E+04, +8.50E+04, +8.30E+04, +8.10E+04, +7.80E+04, +7.50E+04, +7.20E+04, +6.60E+04, +6.00E+04, +5.30E+04, +4.50E+04, +3.70E+04, +2.90E+04, +2.00E+04, +0.00E+00, -3.00E+04, -3.20E+04, -3.30E+04, -3.30E+04, -3.30E+04, -3.20E+04, -3.20E+04, -3.10E+04, -3.10E+04, -3.00E+04, -3.00E+04, -3.00E+04, -2.90E+04, -2.90E+04, -2.80E+04, -2.80E+04, -2.80E+04, -2.80E+04, -2.90E+04, -2.90E+04, -3.00E+04, -3.00E+04, -3.00E+04, -3.10E+04, -3.20E+04, -3.30E+04, -3.50E+04, -3.60E+04, -3.80E+04, -4.00E+04, -4.20E+04, -4.30E+04, -4.40E+04, -4.40E+04, -4.30E+04, -4.10E+04, -3.80E+04, -3.20E+04, -2.10E+04, +0.00E+00};
      if (presetID == 6) ForceLaws[idxMat]= std::vector<float>{1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.00E+08, +0.00E+00, -1.50E+08, -2.30E+08, -1.50E+08, +0.00E+00, +1.10E+08, +1.00E+08, +5.00E+07, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};
    }

    // Compute the effective range of the force law ignoring the zero tail
    constexpr float tailTol= 1.e-2f;
    int tailStartIdx= 0;
    for (int k= 0; k < (int)ForceLaws[idxMat].size(); k++)
      if (std::abs(ForceLaws[idxMat][k]) > std::abs(tailTol * ForceLaws[idxMat][0]))
        tailStartIdx= k;
    ForceLawRanges[idxMat]= ForceLawSteps[idxMat] * (float)(tailStartIdx + 1) * D.UI[LatticePitch____].F();

    // Normalize the force law
    const float BaseVal= ForceLaws[idxMat][0];
    for (int k= 0; k < (int)ForceLaws[idxMat].size(); k++)
      ForceLaws[idxMat][k]/= BaseVal;

    // Scale the force law
    for (int k= 0; k < (int)ForceLaws[idxMat].size(); k++)
      ForceLaws[idxMat][k]*= D.UI[ForceLawScale___].F();
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("ForceLawsT %f\n", Timer::PopTimer());
}


void ParticForceLaw::ComputeBuckets() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Get param
  const int bucketCount= std::max(D.UI[BucketsCount____].I(), 1);
  const int bucketCapacity= std::max(D.UI[BucketsCapacity_].I(), 1);

  // Compute the appropriate voxel size and grid resolution
  const float boxDX= std::max((float)(D.boxMax[0] - D.boxMin[0]), D.UI[LatticePitch____].F());
  const float boxDY= std::max((float)(D.boxMax[1] - D.boxMin[1]), D.UI[LatticePitch____].F());
  const float boxDZ= std::max((float)(D.boxMax[2] - D.boxMin[2]), D.UI[LatticePitch____].F());
  const float stepSize= std::pow(boxDX * boxDY * boxDZ, 1.0f / 3.0f) / std::pow((float)bucketCount, 1.0f / 3.0f);
  nX= std::max(1, (int)std::round(boxDX / stepSize));
  nY= std::max(1, (int)std::round(boxDY / stepSize));
  nZ= std::max(1, (int)std::round(boxDZ / stepSize));

  // Allocate/initialize the spatial partition if needed
  int nXOld, nYOld, nZOld, nBOld;
  Field::GetFieldDimensions(Buckets, nXOld, nYOld, nZOld, nBOld);
  if (nXOld != nX || nYOld != nY || nZOld != nZ || (int)Buckets[0][0][0].capacity() != bucketCapacity)
    Buckets= Field::AllocField4D(nX, nY, nZ, bucketCapacity, -1);
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        Buckets[x][y][z].clear();

  // Fill the spatial partition
  for (int k= 0; k < (int)Pos.size(); k++) {
    int idxX, idxY, idxZ;
    GetBucketIdx(Pos[k], idxX, idxY, idxZ);
    if (idxX >= 0 && idxX < nX && idxY >= 0 && idxY < nY && idxZ >= 0 && idxZ < nZ)
      if ((int)Buckets[idxX][idxY][idxZ].size() < bucketCapacity)
        Buckets[idxX][idxY][idxZ].push_back(k);
  }

  // Plot spatial partition occupancy
  if (D.Plot.size() < 1) D.Plot.resize(1);
  D.Plot[0].name= "Buckets";
  D.Plot[0].val.resize(nX * nY * nZ + 2);
  D.Plot[0].val[nX * nY * nZ + 0]= 0;
  D.Plot[0].val[nX * nY * nZ + 1]= bucketCapacity;
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        D.Plot[0].val[x * nY * nZ + y * nZ + z]= (double)Buckets[x][y][z].size();

  if (D.UI[VerboseLevel____].I() >= 1) printf("SpatialSortT %f\n", Timer::PopTimer());
}


void ParticForceLaw::GetBucketIdx(const Vec::Vec3<float>& iPos, int& oIdxX, int& oIdxY, int& oIdxZ) {
  oIdxX= (iPos[0] == D.boxMax[0]) ? (nX - 1) : ((int)std::floor((float)nX * (iPos[0] - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0])));
  oIdxY= (iPos[1] == D.boxMax[1]) ? (nY - 1) : ((int)std::floor((float)nY * (iPos[1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1])));
  oIdxZ= (iPos[2] == D.boxMax[2]) ? (nZ - 1) : ((int)std::floor((float)nZ * (iPos[2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2])));
}


void ParticForceLaw::ComputeForces() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset forces
  std::fill(For.begin(), For.end(), Vec::Vec3<float>(0.0, 0.0, 0.0));

  // Precompute values
  const float surfArea= 4.0f * std::numbers::pi * D.UI[LatticePitch____].F() * D.UI[LatticePitch____].F();
  float maxForceLawRange= 0.0;
  std::vector<float> ForceLawRangesSqr= ForceLawRanges;
  for (int idxMat= 0; idxMat < (int)ForceLaws.size(); idxMat++) {
    maxForceLawRange= std::max(maxForceLawRange, ForceLawRanges[idxMat]);
    ForceLawRangesSqr[idxMat]= ForceLawRanges[idxMat] * ForceLawRanges[idxMat];
  }

  // Compute the spatial partition for linear neighbor search
  ComputeBuckets();

// Compute particle forces
#pragma omp parallel for
  for (int k0= 0; k0 < (int)Pos.size(); k0++) {
    if (Pos[k0][0] < D.boxMin[0] || Pos[k0][0] > D.boxMax[0] ||
        Pos[k0][1] < D.boxMin[1] || Pos[k0][1] > D.boxMax[1] ||
        Pos[k0][2] < D.boxMin[2] || Pos[k0][2] > D.boxMax[2]) continue;
    // Get range to check in spatial partition
    int idxXBeg, idxYBeg, idxZBeg, idxXEnd, idxYEnd, idxZEnd;
    Vec::Vec3<float> vecOffset(maxForceLawRange, maxForceLawRange, maxForceLawRange);
    GetBucketIdx(Pos[k0] - vecOffset, idxXBeg, idxYBeg, idxZBeg);
    GetBucketIdx(Pos[k0] + vecOffset, idxXEnd, idxYEnd, idxZEnd);
    // Check range in spatial partition
    for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, nX - 1); x++) {
      for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, nY - 1); y++) {
        for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, nZ - 1); z++) {
          // Check candidate particles
          for (int k1 : Buckets[x][y][z]) {
            // Skip if invalid neighbor
            if (k0 == k1) continue;
            if (BCVel[k0] != 0 && BCVel[k1] != 0) continue;
            if (BCPos[k0] != 0 && BCPos[k1] != 0) continue;
            // Precompute distances
            const Vec::Vec3<float> distVec= Pos[k0] - Pos[k1];
            const float distSquared= distVec.normSquared();
            if (distSquared > ForceLawRangesSqr[Mat[k0]] && distSquared > ForceLawRangesSqr[Mat[k1]]) continue;
            const float distVal= std::sqrt(distSquared);
            // Get linear interpolation of force laws for the given distance
            float forceVal= 0.0f;
            std::vector<int> matIndices(1, k0);
            if (k0 != k1) matIndices.push_back(k1);
            for (int k : matIndices) {
              const float valFloat= distVal / (ForceLawSteps[Mat[k]] * D.UI[LatticePitch____].F());
              const int low= std::min(std::max((int)std::floor(valFloat), 0), (int)ForceLaws[Mat[k]].size() - 1);
              const int upp= std::min(std::max(low + 1, 0), (int)ForceLaws[Mat[k]].size() - 1);
              const float ratio= valFloat - (float)low;
              if (k0 == k1) forceVal= ((1.0f - ratio) * ForceLaws[Mat[k]][low] + (ratio)*ForceLaws[Mat[k]][upp]);
              else forceVal+= 0.5f * ((1.0f - ratio) * ForceLaws[Mat[k]][low] + (ratio)*ForceLaws[Mat[k]][upp]);
            }
            // Apply inter-particle force
            For[k0]+= forceVal * surfArea * distVec / distVal;
            // Get radial velocity of particle pair
            const float radialVel= (Vel[k0] - Vel[k1]).dot(distVec / distVal);
            // Apply inter-particle damping proportional to radial velocity
            For[k0]-= (1.0f - distVal / maxForceLawRange) * D.UI[DampingRadRel___].F() * D.UI[ForceLawScale___].F() * radialVel * surfArea * distVec / distVal;
            // TODO find better scaling method ? Current linear decrease up to forcereach means it depends on the furthest reaching material
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
  const float surfArea= 4.0f * std::numbers::pi * D.UI[LatticePitch____].F() * D.UI[LatticePitch____].F();

  // Apply boundary conditions via force controller
  for (int k= 0; k < (int)Pos.size(); k++) {
    if (BCPos[k] != 0) {
      // xt+1 = xt + dt*vt + dt*dt*ft/m
      // xt+1 == xbc  =>  xt + dt*vt + dt*dt*ft/m == xbc
      //              =>  dt*dt*ft/m              == xbc - xt - dt*vt
      //              =>  ft                      == m*(xbc - xt - dt*vt)/dt*dt
      Vec::Vec3<float> ErrVec= Ref[k] - Pos[k];
      For[k]+= D.UI[BCPosCoeff______].F() * surfArea * (ErrVec.normalized() * ErrVec.normSquared() - dt * Vel[k]) / (dt * dt);
    }
    else if (BCVel[k] != 0) {
      // vt+1 = vt + dt*ft/m
      // vt+1 == vbc  =>  vt + dt*ft/m == vbc
      //              =>  dt*ft/m      == vbc - vt
      //              =>  ft           == m*(vbc-vt)/dt
      Vec::Vec3<float> ErrVec= ((BCVel[k] < 0) ? (BCVelVecNega) : (BCVelVecPosi)) - Vel[k];
      For[k]+= D.UI[BCVelCoeff______].F() * surfArea * ErrVec.normalized() * ErrVec.normSquared() / dt;
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
  const float particleMass= D.UI[MaterialDensity_].F() * std::pow(D.UI[LatticePitch____].F(), 3.0f);

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
      Acc[k]= For[k] / particleMass;                           // at = ft / m
      Vel[k]+= dt * Acc[k];                                    // vt+1 = vt + dt * at
      if (!D.UI[UseForceControl_].B() && BCVel[k] != 0)        // Boundary conditions check
        Vel[k]= (BCVel[k] < 0) ? BCVelVecNega : BCVelVecPosi;  // Overwrite velocity
      Pos[k]+= dt * Vel[k];                                    // xt+1 = xt + dt * vt+1
      if (!D.UI[UseForceControl_].B() && BCPos[k] != 0) {      // Boundary conditions check
        Pos[k]= Ref[k];                                        // Overwrite Position
        Vel[k]= Vec::Vec3<float>{0.0f, 0.0f, 0.0f};            // Reset Velocity
      }
      if (D.UI[DampingVelRel___].B())
        Vel[k]= (1.0f - D.UI[DampingVelRel___].F()) * Vel[k];  // Dampen Velocity
    }
  }
  else {
    // Velocity Verlet integration - position update
    for (int k= 0; k < (int)Pos.size(); k++) {
      if (!D.UI[UseForceControl_].B() && BCPos[k] != 0) {  // Boundary conditions check
        Pos[k]= Ref[k];                                    // Overwrite Position
        Vel[k]= Vec::Vec3<float>{0.0f, 0.0f, 0.0f};        // Reset Velocity
      }
      else {
        Pos[k]+= dt * Vel[k] + 0.5f * dt * dt * Acc[k];  // xt+1 = xt + dt * vt + 0.5 * dt * dt * at
      }
      if (D.UI[DampingVelRel___].B())
        Vel[k]= (1.0f - D.UI[DampingVelRel___].F()) * Vel[k];  // Dampen Velocity
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
        Acc[k]= 0.0f;                                          // Reset acceleration
      }
      else {
        const Vec::Vec3<float> oldAcc= Acc[k];   // Store previous acceleration
        Acc[k]= For[k] / particleMass;           // at+1 = ft+1 / m
        Vel[k]+= 0.5f * dt * (oldAcc + Acc[k]);  // vt+1 = vt + 0.5 * dt * (at + at+1)
      }
    }
  }

  // Advance time
  SimTime+= dt;

  // Set particle color with gradual decay
  for (int k= 0; k < (int)Pos.size(); k++) {
    float r= 0.5f, g= 0.5f, b= 0.5f;
    if (D.UI[ColorMode_______].I() == 0) {
      if (Mat[k] == 0) r= 1.0f;
      if (Mat[k] == 1) g= 1.0f;
      if (Mat[k] == 2) b= 1.0f;
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
    if (D.UI[ColorMode_______].I() == 3) Colormap::RatioToJetBrightSmooth(For[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
    Col[k]= (1.0f - D.UI[ColorDecay______].F()) * Col[k] + D.UI[ColorDecay______].F() * Vec::Vec3<float>(r, g, b);
  }

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

void ParticForceLaw::ComputeMetaballs() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  MetaballIsUpdated= true;
  Verts.clear();
  Tris.clear();

  if (D.UI[MetaballVoxSize_].F() <= 0.0f) return;

  const float metaballSize= D.UI[LatticePitch____].F();
  const int tmpNX= std::max((int)std::round((D.boxMax[0] - D.boxMin[0]) / D.UI[MetaballVoxSize_].D()), 1);
  const int tmpNY= std::max((int)std::round((D.boxMax[1] - D.boxMin[1]) / D.UI[MetaballVoxSize_].D()), 1);
  const int tmpNZ= std::max((int)std::round((D.boxMax[2] - D.boxMin[2]) / D.UI[MetaballVoxSize_].D()), 1);
  std::vector<std::vector<std::vector<double>>> tmpField= Field::AllocField3D(tmpNX, tmpNY, tmpNZ, 0.0);
  double stepX, stepY, stepZ;
  BoxGrid::GetVoxelSizes(tmpNX, tmpNY, tmpNZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ);
  for (int k= 0; k < (int)Pos.size(); k++) {
    const int idxXBeg= (int)std::floor((float)tmpNX * (Pos[k][0] - 2.0 * metaballSize - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxYBeg= (int)std::floor((float)tmpNY * (Pos[k][1] - 2.0 * metaballSize - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZBeg= (int)std::floor((float)tmpNZ * (Pos[k][2] - 2.0 * metaballSize - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    const int idxXEnd= (int)std::floor((float)tmpNX * (Pos[k][0] + 2.0 * metaballSize - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxYEnd= (int)std::floor((float)tmpNY * (Pos[k][1] + 2.0 * metaballSize - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZEnd= (int)std::floor((float)tmpNZ * (Pos[k][2] + 2.0 * metaballSize - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, tmpNX - 1); x++) {
      for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, tmpNY - 1); y++) {
        for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, tmpNZ - 1); z++) {
          const Vec::Vec3<float> pos(D.boxMin[0] + (x + 0.5f) * stepX, D.boxMin[1] + (y + 0.5f) * stepY, D.boxMin[2] + (z + 0.5f) * stepZ);
          if ((pos - Pos[k]).normSquared() < (2.0 * metaballSize) * (2.0 * metaballSize)) {
            tmpField[x][y][z]+= std::max(0.0, 1.0 - (pos - Pos[k]).norm() / metaballSize);
          }
        }
      }
    }
  }

  MarchingCubes::ComputeMarchingCubes(D.UI[MetaballIsoval__].F(), D.boxMin, D.boxMax, tmpField, Verts, Tris);

  if (D.UI[VerboseLevel____].I() >= 1) printf("MetaballT %f\n", Timer::PopTimer());
}
