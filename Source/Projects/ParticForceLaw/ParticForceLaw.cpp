#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <format>
#include <numbers>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

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
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("DomainX_________", 1.0));    // Simulation domain size
    D.UI.push_back(ParamUI("DomainY_________", 1.0));    // Simulation domain size
    D.UI.push_back(ParamUI("DomainZ_________", 1.0));    // Simulation domain size
    D.UI.push_back(ParamUI("ScenarioPreset__", 4));      // Scenario choice, 0 load file, 1+ hard coded scenarios
    D.UI.push_back(ParamUI("Scenario2DID____", 0));      // BMP file to load
    D.UI.push_back(ParamUI("Scenario2DThick_", 0.5));    // Relative thickness of loaded file wrt domain
    D.UI.push_back(ParamUI("LatticePitch____", 0.04));   // Charateristic lengthscale of particle cloud, controls the number of particles
    D.UI.push_back(ParamUI("LatticePattern__", 2));      // Choose the particle cloud pattern
    D.UI.push_back(ParamUI("ConstrainDim2D__", 0));      // Optionally constrain movement to plane orthogonal to X=1, Y=2, X=3 axis
    D.UI.push_back(ParamUI("______________00", NAN));    //
    D.UI.push_back(ParamUI("StepsPerDraw____", 1));      // Number of simulation steps between each draw call
    D.UI.push_back(ParamUI("TimeStep________", 1.e-3));  // Physical timestep of explicit integration
    D.UI.push_back(ParamUI("BucketCapacity__", 40));     // Chosen bucket capacity
    D.UI.push_back(ParamUI("BucketFillCoeff_", 4.0));    // Scaling parameter of spatial partition
    D.UI.push_back(ParamUI("DampingRadRel___", 0.1));    // Physical interparticle damping proportional to difference of radial velocities
    D.UI.push_back(ParamUI("DampingVelRel___", 0.0));    // Unphysical absolute velocity damping for stability
    D.UI.push_back(ParamUI("MaterialDensity_", 200.0));  // Density of simulated material
    D.UI.push_back(ParamUI("______________01", NAN));    //
    D.UI.push_back(ParamUI("BCVelX__________", 0.0));    // Velocity for particles with associated BC
    D.UI.push_back(ParamUI("BCVelY__________", 0.0));    // Velocity for particles with associated BC
    D.UI.push_back(ParamUI("BCVelZ__________", 1.0));    // Velocity for particles with associated BC
    D.UI.push_back(ParamUI("BCForX__________", 0.0));    // External force for particles with associated BC
    D.UI.push_back(ParamUI("BCForY__________", 0.0));    // External force for particles with associated BC
    D.UI.push_back(ParamUI("BCForZ__________", 1.0));    // External force for particles with associated BC
    D.UI.push_back(ParamUI("UseForceControl_", 0));      // Enable force controller instead of simple overwrite
    D.UI.push_back(ParamUI("BCPosCoeff______", 1.0));    // Scaling coeff for force controller on position BC
    D.UI.push_back(ParamUI("BCVelCoeff______", 1.0));    // Scaling coeff for force controller on velocity BC
    D.UI.push_back(ParamUI("______________02", NAN));    //
    D.UI.push_back(ParamUI("ForceLawPresetA_", 0));      // Chosen force law for material 1
    D.UI.push_back(ParamUI("ForceLawPresetB_", 1));      // Chosen force law for material 2
    D.UI.push_back(ParamUI("ForceLawScaleA__", 100.0));  // Force law scaling for material 1
    D.UI.push_back(ParamUI("ForceLawScaleB__", 100.0));  // Force law scaling for material 2
    D.UI.push_back(ParamUI("______________03", NAN));    //
    D.UI.push_back(ParamUI("ForceLawA_0_00__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_0_80__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_0_90__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_0_95__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_00__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_05__", -1.0));   // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_10__", -1.0));   // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_20__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_30__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_40__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_50__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_2_00__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_2_50__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_3_00__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("______________04", NAN));    //
    D.UI.push_back(ParamUI("ForceLawB_0_00__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_0_80__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_0_90__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_0_95__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_00__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_05__", -1.0));   // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_10__", -1.0));   // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_20__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_30__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_40__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_50__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_2_00__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_2_50__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_3_00__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("______________05", NAN));    //
    D.UI.push_back(ParamUI("MetaballVoxSize_", 0.04));   // Voxel size for Metaball and Marching Cubes isosurface
    D.UI.push_back(ParamUI("MetaballIsoval__", 0.5));    // Isovalue for Marching Cubes
    D.UI.push_back(ParamUI("ColorMode_______", 3));      // Show material, boundary conditions, velocity, force magnitude, velocity magnitude
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));    // Color factor for modes supporting it
    D.UI.push_back(ParamUI("VisuScale_______", 0.5));    // Size scaling of display elements
    D.UI.push_back(ParamUI("VisuSimple______", 0));      // Toggle for simplified draw mode
    D.UI.push_back(ParamUI("VisuHideOOB_____", 0));      // Hide particles out of bounds
    D.UI.push_back(ParamUI("VisuMinNeighbor_", 0));      // Hide particles with too few neighbors
    D.UI.push_back(ParamUI("______________06", NAN));    //
    D.UI.push_back(ParamUI("TestParamMIP_0__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_1__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_2__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_3__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_4__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_5__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("VerboseLevel____", 0));      // Verbose mode
  }

  if (D.UI.size() != VerboseLevel____ + 1) {
    printf("[ERROR] Invalid parameter count in UI\n");
  }

  RunID= 0;

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
  if (D.UI[ForceLawPresetA_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawPresetB_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawScaleA__].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawScaleB__].hasChanged()) isRefreshed= false;
  for (int idxParam= ForceLawA_0_00__; idxParam <= ForceLawA_3_00__; idxParam++)
    if (D.UI[idxParam].hasChanged()) isRefreshed= false;
  for (int idxParam= ForceLawB_0_00__; idxParam <= ForceLawB_3_00__; idxParam++)
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
  For.clear();
  Col.clear();
  Mat.clear();
  ForceMag.clear();
  Neighbors.clear();
  Sensor.clear();
  BCPos.clear();
  BCVel.clear();
  BCFor.clear();
  MetaballIsUpdated= false;
  Verts.clear();
  Tris.clear();

  // Get domain dimensions
  D.boxMin= {0.5 - 0.5 * std::max(D.UI[DomainX_________].D(), D.UI[LatticePitch____].D()),
             0.5 - 0.5 * std::max(D.UI[DomainY_________].D(), D.UI[LatticePitch____].D()),
             0.5 - 0.5 * std::max(D.UI[DomainZ_________].D(), D.UI[LatticePitch____].D())};
  D.boxMax= {0.5 + 0.5 * std::max(D.UI[DomainX_________].D(), D.UI[LatticePitch____].D()),
             0.5 + 0.5 * std::max(D.UI[DomainY_________].D(), D.UI[LatticePitch____].D()),
             0.5 + 0.5 * std::max(D.UI[DomainZ_________].D(), D.UI[LatticePitch____].D())};

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
  for (int idxMat= 0; idxMat < (int)ForceLaw.size(); idxMat++) {
    D.Plot[3 + idxMat].name= "ForceLaw";
    D.Plot[3 + idxMat].isSymmetric= true;
    if (idxMat > 0) D.Plot[3 + idxMat].isSameRange= true;
    D.Plot[3 + idxMat].val.resize(ForceLaw[idxMat].size());
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      D.Plot[3 + idxMat].val[k]= ForceLaw[idxMat][k];
  }
}


// Handle keypress
void ParticForceLaw::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Handle mouse action
void ParticForceLaw::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
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
  D.Plot[2].name= "PE";
  if (D.Plot[1].val.size() < 5000) {
    D.Plot[1].val.reserve(5000);
    D.Plot[2].val.reserve(5000);
    D.Plot[1].val.push_back(0.0f);
    D.Plot[2].val.push_back(0.0f);
    for (int k= 0; k < (int)Pos.size(); k++) {
      D.Plot[1].val[D.Plot[1].val.size() - 1]+= Vel[k].normSquared();
      D.Plot[2].val[D.Plot[2].val.size() - 1]+= ForceMag[k];
    }
  }
}


// Draw the project
void ParticForceLaw::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Set particle color
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
    if (D.UI[ColorMode_______].I() == 3) {
      Colormap::RatioToJetBrightSmooth(ForceMag[k] * D.UI[ColorFactor_____].F(), r, g, b);
    }
    if (D.UI[ColorMode_______].I() == 4) {
      Colormap::RatioToBlackBody(Vel[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
    }
    Col[k].set(r, g, b);
  }

  // Display particles
  if (D.displayMode1) {
    if (D.UI[VisuSimple______].B()) {
      glPointSize(1000.0f * D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F());
      glBegin(GL_POINTS);
      for (int k= 0; k < (int)Pos.size(); k++) {
        if (D.UI[VisuHideOOB_____].B() &&
            (Pos[k][0] < D.boxMin[0] || Pos[k][0] > D.boxMax[0] ||
             Pos[k][1] < D.boxMin[1] || Pos[k][1] > D.boxMax[1] ||
             Pos[k][2] < D.boxMin[2] || Pos[k][2] > D.boxMax[2])) continue;
        if (Neighbors[k] < D.UI[VisuMinNeighbor_].I()) continue;
        glColor3fv(Col[k].array());
        glVertex3fv(Pos[k].array());
      }
      glEnd();
    }
    else {
      glEnable(GL_LIGHTING);
      for (int k= 0; k < (int)Pos.size(); k++) {
        if (D.UI[VisuHideOOB_____].B() &&
            (Pos[k][0] < D.boxMin[0] || Pos[k][0] > D.boxMax[0] ||
             Pos[k][1] < D.boxMin[1] || Pos[k][1] > D.boxMax[1] ||
             Pos[k][2] < D.boxMin[2] || Pos[k][2] > D.boxMax[2])) continue;
        if (Neighbors[k] < D.UI[VisuMinNeighbor_].I()) continue;
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
    BoxGrid::GetVoxelSizes(Buckets.nX, Buckets.nY, Buckets.nZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ);
    const int bucketCapacity= std::max(D.UI[BucketCapacity__].I(), 1);
    // Set transformation
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5f * (float)stepX, D.boxMin[1] + 0.5f * (float)stepY, D.boxMin[2] + 0.5f * (float)stepZ);
    glScalef((float)stepX, (float)stepY, (float)stepZ);
    for (int x= 0; x < Buckets.nX; x++) {
      for (int y= 0; y < Buckets.nY; y++) {
        for (int z= 0; z < Buckets.nZ; z++) {
          // Color by occupancy
          float r= 0.5f, g= 0.5f, b= 0.5f;
          Colormap::RatioToJetBrightSmooth((float)Buckets.at(x, y, z).size() / (float)bucketCapacity, r, g, b);
          glColor3f(r, g, b);
          // Draw wire box
          glPushMatrix();
          glTranslatef((float)x, (float)y, (float)z);
          glutWireCube(0.95);
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
  D.Status.clear();
  D.Status.resize(4);
  D.Status[0]= std::format("NbBuckets:{}", Buckets.nXYZ);
  D.Status[1]= std::format("NbParticles:{}", (int)Pos.size());
  D.Status[2]= std::format("SimTime:{:.6f}ms", SimTime);
  if (BucketOverflown) D.Status[3]= std::string{"BUCKET OVERFLOW"};

  if (D.UI[VerboseLevel____].I() >= 1) printf("DrawT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the base point cloud
  oPointCloud.clear();

  // Check parameters
  if (D.UI[LatticePitch____].F() <= 0.0f) return;

  // Regular cubic lattice patterns
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
  // HCP pattern with layers along X
  else if (D.UI[LatticePattern__].I() == 3) {
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * std::sqrt(6.0f) / 3.0f));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F()));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * 0.5f * std::sqrt(3.0f)));
    for (int x= 0; x < cellNbX + 1; x++)
      for (int y= 0; y < cellNbY; y++)
        for (int z= 0; z < cellNbZ; z++)
          oPointCloud.push_back(0.5f * D.UI[LatticePitch____].F() *
                                Vec::Vec3<float>(float(x) * 2.0f * std::sqrt(6.0f) / 3.0f,
                                                 2.0f * float(y) + float((z + x) % 2),
                                                 std::sqrt(3.0f) * (float(z) + float(x % 2) / 3.0f)));
  }
  // HCP pattern with layers along Y
  else if (D.UI[LatticePattern__].I() == 4) {
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * 0.5f * std::sqrt(3.0f)));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * std::sqrt(6.0f) / 3.0f));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F()));
    for (int x= 0; x < cellNbX; x++)
      for (int y= 0; y < cellNbY + 1; y++)
        for (int z= 0; z < cellNbZ; z++)
          oPointCloud.push_back(0.5f * D.UI[LatticePitch____].F() *
                                Vec::Vec3<float>(std::sqrt(3.0f) * (float(x) + float(y % 2) / 3.0f),
                                                 float(y) * 2.0f * std::sqrt(6.0f) / 3.0f,
                                                 2.0f * float(z) + float((x + y) % 2)));
  }
  // HCP pattern with layers along Z
  else if (D.UI[LatticePattern__].I() == 5) {
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F()));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * 0.5f * std::sqrt(3.0f)));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * std::sqrt(6.0f) / 3.0f));
    for (int x= 0; x < cellNbX; x++)
      for (int y= 0; y < cellNbY; y++)
        for (int z= 0; z < cellNbZ + 1; z++)
          oPointCloud.push_back(0.5f * D.UI[LatticePitch____].F() *
                                Vec::Vec3<float>(2.0f * float(x) + float((y + z) % 2),
                                                 std::sqrt(3.0f) * (float(y) + float(z % 2) / 3.0f),
                                                 float(z) * 2.0f * std::sqrt(6.0f) / 3.0f));
  }
  // Poisson sphere sampling
  else if (D.UI[LatticePattern__].I() == 6) {
    ComputeBuckets();
    const int bucketCapacity= std::max(D.UI[BucketCapacity__].I(), 1);
    int failStreak= 0;
    const int maxAttempts= 1000;
    const float relMinDist= 0.9f;
    while (failStreak < maxAttempts) {
      Vec::Vec3<float> candidate(Random::Val(D.boxMin[0], D.boxMax[0]), Random::Val(D.boxMin[1], D.boxMax[1]), Random::Val(D.boxMin[2], D.boxMax[2]));
      if (D.UI[ConstrainDim2D__].I() == 1) candidate[0]= D.boxMin[0] + 0.5 * (D.boxMax[0] - D.boxMin[0]);
      if (D.UI[ConstrainDim2D__].I() == 2) candidate[1]= D.boxMin[1] + 0.5 * (D.boxMax[1] - D.boxMin[1]);
      if (D.UI[ConstrainDim2D__].I() == 3) candidate[2]= D.boxMin[2] + 0.5 * (D.boxMax[2] - D.boxMin[2]);
      bool keep= true;
      // Get range to check in spatial partition
      int idxXBeg, idxYBeg, idxZBeg, idxXEnd, idxYEnd, idxZEnd;
      Vec::Vec3<float> vecOffset(D.UI[LatticePitch____].F(), D.UI[LatticePitch____].F(), D.UI[LatticePitch____].F());
      GetBucketIdx(candidate - vecOffset, idxXBeg, idxYBeg, idxZBeg);
      GetBucketIdx(candidate + vecOffset, idxXEnd, idxYEnd, idxZEnd);
      // Check range in spatial partition
      for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, Buckets.nX - 1) && keep; x++) {
        for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, Buckets.nY - 1) && keep; y++) {
          for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, Buckets.nZ - 1) && keep; z++) {
            // Check candidate particles
            for (int k : Buckets.at(x, y, z)) {
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
        if (idxX >= 0 && idxX < Buckets.nX && idxY >= 0 && idxY < Buckets.nY && idxZ >= 0 && idxZ < Buckets.nZ)
          if ((int)Buckets.at(idxX, idxY, idxZ).size() < bucketCapacity)
            Buckets.at(idxX, idxY, idxZ).push_back((int)oPointCloud.size() - 1);
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
    if (D.UI[Scenario2DID____].I() == 3) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/Auxetic.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 4) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/Logo.bmp", imageRGBA, false);
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
    // Ball blasting through Bimaterial cylinder
    else if (D.UI[ScenarioPreset__].I() == 11) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.15f)).norm() < 0.10f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Mat.push_back(0);
      }
      else if (RelPos[1] > 0.1f && RelPos[1] < 0.9f && (RelPos - Vec::Vec3<float>(0.5f, RelPos[1], 0.5f)).norm() < 0.05f) {
        Pos.push_back(iPointCloud[k]);
        Mat.push_back(0);
      }
      else if (RelPos[1] > 0.1f && RelPos[1] < 0.9f && (RelPos - Vec::Vec3<float>(0.5f, RelPos[1], 0.5f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Mat.push_back(1);
      }
    }

    // Fill the missing default values
    if (Mat.size() < Pos.size()) Mat.push_back(0);
    if (ForceMag.size() < Pos.size()) ForceMag.push_back(0.0);
    if (Neighbors.size() < Pos.size()) Neighbors.push_back(1);
    if (Sensor.size() < Pos.size()) Sensor.push_back(0);
    if (BCPos.size() < Pos.size()) BCPos.push_back(0);
    if (BCVel.size() < Pos.size()) BCVel.push_back(0);
    if (BCFor.size() < Pos.size()) BCFor.push_back(0);

    if (Ref.size() < Pos.size()) Ref.push_back(Pos[Pos.size() - 1]);
    if (Vel.size() < Pos.size()) Vel.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (For.size() < Pos.size()) For.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Col.size() < Pos.size()) Col.push_back(Vec::Vec3<float>(0.5f, 0.5f, 0.5f));
  }

  SimTime= 0.0f;
  RunID++;

  if (D.UI[VerboseLevel____].I() >= 1) printf("ScenarioT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildForceLaws() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the force laws
  constexpr int nbForceLaws= 2;
  ForceLaw.clear();
  ForceLaw.resize(nbForceLaws);
  ForceLawStep.clear();
  ForceLawStep.resize(nbForceLaws);
  ForceLawRange.clear();
  ForceLawRange.resize(nbForceLaws);

  // Create the force laws
  for (int idxMat= 0; idxMat < (int)ForceLaw.size(); idxMat++) {
    int presetID= -1;
    float forceLawScale= 1.0f;
    if (idxMat == 0) {
      presetID= D.UI[ForceLawPresetA_].I();
      forceLawScale= D.UI[ForceLawScaleA__].F();
    }
    if (idxMat == 1) {
      presetID= D.UI[ForceLawPresetB_].I();
      forceLawScale= D.UI[ForceLawScaleB__].F();
    }
    // Custom force law
    if (presetID == 0) {
      BuildForceLawPolyline(D.UI[ForceLawA_0_00__].D(), D.UI[ForceLawA_0_80__].D(), D.UI[ForceLawA_0_90__].D(), D.UI[ForceLawA_0_95__].D(),
                            D.UI[ForceLawA_1_00__].D(), D.UI[ForceLawA_1_05__].D(), D.UI[ForceLawA_1_10__].D(), D.UI[ForceLawA_1_20__].D(),
                            D.UI[ForceLawA_1_30__].D(), D.UI[ForceLawA_1_40__].D(), D.UI[ForceLawA_1_50__].D(), D.UI[ForceLawA_2_00__].D(),
                            D.UI[ForceLawA_2_50__].D(), D.UI[ForceLawA_3_00__].D(), idxMat);
    }
    else if (presetID == 1) {
      BuildForceLawPolyline(D.UI[ForceLawB_0_00__].D(), D.UI[ForceLawB_0_80__].D(), D.UI[ForceLawB_0_90__].D(), D.UI[ForceLawB_0_95__].D(),
                            D.UI[ForceLawB_1_00__].D(), D.UI[ForceLawB_1_05__].D(), D.UI[ForceLawB_1_10__].D(), D.UI[ForceLawB_1_20__].D(),
                            D.UI[ForceLawB_1_30__].D(), D.UI[ForceLawB_1_40__].D(), D.UI[ForceLawB_1_50__].D(), D.UI[ForceLawB_2_00__].D(),
                            D.UI[ForceLawB_2_50__].D(), D.UI[ForceLawB_3_00__].D(), idxMat);
    }
    // Hard coded force laws from MIT paper
    else {
      ForceLawStep[idxMat]= 0.05f;
      ForceLaw[idxMat]= std::vector<float>{1.0f};
      //                                                      0.00,      0.05,      0.10,      0.15,      0.20,      0.25,      0.30,      0.35,      0.40,      0.45,      0.50,      0.55,      0.60,      0.65,      0.70,      0.75,      0.80,      0.85,      0.90,      0.95,      1.00,      1.05,      1.10,      1.15,      1.20,      1.25,      1.30,      1.35,      1.40,      1.45,      1.50,      1.55,      1.60,      1.65,      1.70,      1.75,      1.80,      1.85,      1.90,      1.95,      2.00,      2.05,      2.10,      2.15,      2.20,      2.25,      2.30,      2.35,      2.40,      2.45,      2.50,      2.55,      2.60,      2.65,      2.70,      2.75,      2.80,      2.85,      2.90,      2.95,      3.00
      if (presetID == 2) ForceLaw[idxMat]= std::vector<float>{1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +5.00E+05, +0.00E+00, -5.00E+05, -1.00E+06, -5.00E+05, +0.00E+00, +5.00E+05, +1.00E+06, +5.00E+05, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Elastic material
      if (presetID == 3) ForceLaw[idxMat]= std::vector<float>{1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +5.00E+05, +0.00E+00, -1.00E+05, +0.00E+00, +7.00E+04, +1.00E+05, +1.20E+05, +8.00E+04, +3.00E+04, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Brittle Material
      if (presetID == 4) ForceLaw[idxMat]= std::vector<float>{1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +9.50E+01, +8.80E+01, +7.00E+01, +5.00E+01, +2.50E+01, +0.00E+00, -2.00E+01, -2.80E+01, -3.00E+01, -2.80E+01, -2.50E+01, -2.00E+01, -1.20E+01, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Viscous material
      if (presetID == 5) ForceLaw[idxMat]= std::vector<float>{3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +2.00E+10, +0.00E+00, -2.00E+10, -2.50E+10, -2.00E+10, +0.00E+00, +2.00E+10, +2.50E+10, +2.00E+10, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Steel AISI 4340
      if (presetID == 6) ForceLaw[idxMat]= std::vector<float>{1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.00E+08, +0.00E+00, -1.50E+08, -2.30E+08, -1.50E+08, +0.00E+00, +1.10E+08, +1.00E+08, +5.00E+07, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Sample force law
      if (presetID == 7) ForceLaw[idxMat]= std::vector<float>{8.80E+04, +8.80E+04, +8.80E+04, +8.70E+04, +8.70E+04, +8.60E+04, +8.50E+04, +8.50E+04, +8.30E+04, +8.10E+04, +7.80E+04, +7.50E+04, +7.20E+04, +6.60E+04, +6.00E+04, +5.30E+04, +4.50E+04, +3.70E+04, +2.90E+04, +2.00E+04, +0.00E+00, -3.00E+04, -3.20E+04, -3.30E+04, -3.30E+04, -3.30E+04, -3.20E+04, -3.20E+04, -3.10E+04, -3.10E+04, -3.00E+04, -3.00E+04, -3.00E+04, -2.90E+04, -2.90E+04, -2.80E+04, -2.80E+04, -2.80E+04, -2.80E+04, -2.90E+04, -2.90E+04, -3.00E+04, -3.00E+04, -3.00E+04, -3.10E+04, -3.20E+04, -3.30E+04, -3.50E+04, -3.60E+04, -3.80E+04, -4.00E+04, -4.20E+04, -4.30E+04, -4.40E+04, -4.40E+04, -4.30E+04, -4.10E+04, -3.80E+04, -3.20E+04, -2.10E+04, +0.00E+00};  // Delrin force law optimized on force displacement curve
    }

    // Compute the effective range of the force law ignoring the zero tail
    int tailStartIdx= 0;
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      if (std::abs(ForceLaw[idxMat][k]) > 0.0f)
        tailStartIdx= k;
    ForceLawRange[idxMat]= ForceLawStep[idxMat] * (float)(tailStartIdx + 1) * D.UI[LatticePitch____].F();

    // Normalize the force law by the first value
    const float BaseVal= ForceLaw[idxMat][0];
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      ForceLaw[idxMat][k]/= BaseVal;

    // Scale the force law by the UI coeff
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      ForceLaw[idxMat][k]*= forceLawScale;
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("ForceLawsT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildForceLawPolyline(const double v0_00, const double v0_80, const double v0_90, const double v0_95,
                                           const double v1_00, const double v1_05, const double v1_10, const double v1_20,
                                           const double v1_30, const double v1_40, const double v1_50, const double v2_00,
                                           const double v2_50, const double v3_00, const int iIdxMat) {
  // Create the force law as a smoothed polyline
  std::vector<std::array<double, 3>> PolylineA;
  PolylineA.push_back(std::array<double, 3>{0.0, 0.00, v0_00});
  PolylineA.push_back(std::array<double, 3>{0.0, 0.80, v0_80});
  PolylineA.push_back(std::array<double, 3>{0.0, 0.90, v0_90});
  PolylineA.push_back(std::array<double, 3>{0.0, 0.95, v0_95});
  PolylineA.push_back(std::array<double, 3>{0.0, 1.00, v1_00});
  std::vector<std::array<double, 3>> PolylineB;
  PolylineB.push_back(std::array<double, 3>{0.0, 1.05, v1_05});
  PolylineB.push_back(std::array<double, 3>{0.0, 1.10, v1_10});
  PolylineB.push_back(std::array<double, 3>{0.0, 1.20, v1_20});
  PolylineB.push_back(std::array<double, 3>{0.0, 1.30, v1_30});
  PolylineB.push_back(std::array<double, 3>{0.0, std::sqrt(2), v1_40});
  if (v1_40 != 0.0 || v1_50 != 0.0 || v2_00 != 0.0 || v2_50 != 0.0 || v3_00 != 0.0) {
    PolylineB.push_back(std::array<double, 3>{0.0, 1.50, v1_50});
    PolylineB.push_back(std::array<double, 3>{0.0, 2.00, v2_00});
    PolylineB.push_back(std::array<double, 3>{0.0, 2.50, v2_50});
    PolylineB.push_back(std::array<double, 3>{0.0, 3.00, v3_00});
  }
  Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineA);
  Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineB);
  std::vector<std::array<double, 3>> Polyline;
  Polyline.insert(Polyline.end(), PolylineA.begin(), PolylineA.end());
  Polyline.insert(Polyline.end(), PolylineB.begin(), PolylineB.end());

  // Sample the force value at fixed distance intervals
  constexpr int nbSamples= 211;  // 3*70+1 chosen to have a sample point at 1.0 and another very close to sqrt(2)
  constexpr float maxReach= 3.0f;
  ForceLaw[iIdxMat].resize(nbSamples);
  ForceLawStep[iIdxMat]= maxReach / float(nbSamples - 1);
  for (int k= 0; k < nbSamples; k++) {
    float dist= ForceLawStep[iIdxMat] * float(k);
    int idxLow= 0, idxUpp= 0;
    for (int idxVert= 0; idxVert < (int)Polyline.size() - 1; idxVert++) {
      idxLow= idxVert;
      idxUpp= std::min(idxVert + 1, (int)Polyline.size() - 1);
      if (Polyline[idxLow][1] <= dist && Polyline[idxUpp][1] >= dist) break;
    }
    double ratio= 0.0;
    if (Polyline[idxUpp][1] - Polyline[idxLow][1] > 0.0)
      ratio= (dist - Polyline[idxLow][1]) / (Polyline[idxUpp][1] - Polyline[idxLow][1]);
    ratio= std::min(std::max(ratio, 0.0), 1.0);
    ForceLaw[iIdxMat][k]= (1.0 - ratio) * Polyline[idxLow][2] + ratio * Polyline[idxUpp][2];
  }
}


void ParticForceLaw::ComputeBuckets() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Compute appropriate bucket grid resolution
  const int bucketCapacity= std::max(D.UI[BucketCapacity__].I(), 1);
  const float particleDensity= D.UI[BucketFillCoeff_].F() / std::pow(D.UI[LatticePitch____].F(), 3.0f);
  const float boxDX= (float)(D.boxMax[0] - D.boxMin[0]);
  const float boxDY= (float)(D.boxMax[1] - D.boxMin[1]);
  const float boxDZ= (float)(D.boxMax[2] - D.boxMin[2]);
  const float theoryBucketCount= particleDensity * boxDX * boxDY * boxDZ / (float)bucketCapacity;
  const float stepSize= std::pow(boxDX * boxDY * boxDZ, 1.0f / 3.0f) / std::pow(theoryBucketCount, 1.0f / 3.0f);
  int const nX= std::max(1, (int)std::round(boxDX / stepSize));
  int const nY= std::max(1, (int)std::round(boxDY / stepSize));
  int const nZ= std::max(1, (int)std::round(boxDZ / stepSize));

  // Allocate and initialize the spatial partition if neede
  if (Buckets.nX != nX || Buckets.nY != nY || Buckets.nZ != nZ || (int)Buckets.at(0, 0, 0).capacity() != bucketCapacity)
    Buckets= Field::Field3<std::vector<int>>(nX, nY, nZ, std::vector<int>(bucketCapacity, -1));
  for (int xyz= 0; xyz < Buckets.nXYZ; xyz++)
    Buckets.at(xyz).clear();
  BucketOverflown= false;

  // Fill the spatial partition
  for (int k= 0; k < (int)Pos.size(); k++) {
    int idxX, idxY, idxZ;
    GetBucketIdx(Pos[k], idxX, idxY, idxZ);
    if (idxX >= 0 && idxX < nX && idxY >= 0 && idxY < nY && idxZ >= 0 && idxZ < nZ) {
      if ((int)Buckets.at(idxX, idxY, idxZ).size() < bucketCapacity)
        Buckets.at(idxX, idxY, idxZ).push_back(k);
      else
        BucketOverflown= true;
    }
  }

  // Plot spatial partition occupancy
  if (D.Plot.size() < 1) D.Plot.resize(1);
  D.Plot[0].name= "Buckets";
  D.Plot[0].val.resize(Buckets.nXYZ + 2);
  D.Plot[0].val[Buckets.nXYZ + 0]= 0;
  D.Plot[0].val[Buckets.nXYZ + 1]= bucketCapacity;
  for (int xyz= 0; xyz < Buckets.nXYZ; xyz++)
    D.Plot[0].val[xyz]= (double)Buckets.at(xyz).size();

  if (D.UI[VerboseLevel____].I() >= 1) printf("SpatialSortT %f\n", Timer::PopTimer());
}


void ParticForceLaw::GetBucketIdx(const Vec::Vec3<float>& iPos, int& oIdxX, int& oIdxY, int& oIdxZ) {
  oIdxX= (iPos[0] == D.boxMax[0]) ? (Buckets.nX - 1) : ((int)std::floor((float)Buckets.nX * (iPos[0] - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0])));
  oIdxY= (iPos[1] == D.boxMax[1]) ? (Buckets.nY - 1) : ((int)std::floor((float)Buckets.nY * (iPos[1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1])));
  oIdxZ= (iPos[2] == D.boxMax[2]) ? (Buckets.nZ - 1) : ((int)std::floor((float)Buckets.nZ * (iPos[2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2])));
}


void ParticForceLaw::ComputeForces() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset vectors
  std::fill(For.begin(), For.end(), Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  std::fill(ForceMag.begin(), ForceMag.end(), 0.0f);
  std::fill(Neighbors.begin(), Neighbors.end(), 0);

  // Precompute values
  const float surfArea= 4.0f * std::numbers::pi * D.UI[LatticePitch____].F() * D.UI[LatticePitch____].F();
  float maxForceLawRange= 0.0;
  std::vector<float> ForceLawRangesSqr= ForceLawRange;
  for (int idxMat= 0; idxMat < (int)ForceLaw.size(); idxMat++) {
    maxForceLawRange= std::max(maxForceLawRange, ForceLawRange[idxMat]);
    ForceLawRangesSqr[idxMat]= ForceLawRange[idxMat] * ForceLawRange[idxMat];
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
    for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, Buckets.nX - 1); x++) {
      for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, Buckets.nY - 1); y++) {
        for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, Buckets.nZ - 1); z++) {
          // Check candidate particles
          for (int k1 : Buckets.at(x, y, z)) {
            // Skip if invalid neighbor
            // if (BCVel[k0] != 0 && BCVel[k1] != 0) continue;
            // if (BCPos[k0] != 0 && BCPos[k1] != 0) continue;
            if (k0 == k1) continue;
            // Precompute distances
            const Vec::Vec3<float> distVec= Pos[k0] - Pos[k1];
            const float distSquared= distVec.normSquared();
            if (distSquared > ForceLawRangesSqr[Mat[k0]] && distSquared > ForceLawRangesSqr[Mat[k1]]) continue;
            const float distVal= std::sqrt(distSquared);
            const Vec::Vec3<float> distVecUnit= distVec / distVal;
            // Get linear interpolation of force laws for the given distance
            float forceVal= 0.0f;
            std::vector<int> matIndices(1, k0);
            if (k0 != k1) matIndices.push_back(k1);
            for (int k : matIndices) {
              const float valFloat= distVal / (ForceLawStep[Mat[k]] * D.UI[LatticePitch____].F());
              const int low= std::min(std::max((int)std::floor(valFloat), 0), (int)ForceLaw[Mat[k]].size() - 1);
              const int upp= std::min(std::max(low + 1, 0), (int)ForceLaw[Mat[k]].size() - 1);
              const float ratio= valFloat - (float)low;
              if (k0 == k1) forceVal= ((1.0f - ratio) * ForceLaw[Mat[k]][low] + (ratio)*ForceLaw[Mat[k]][upp]);
              else forceVal+= 0.5f * ((1.0f - ratio) * ForceLaw[Mat[k]][low] + (ratio)*ForceLaw[Mat[k]][upp]);
            }
            // Apply inter-particle force
            ForceMag[k0]+= std::abs(forceVal) * surfArea;
            For[k0]+= forceVal * surfArea * distVecUnit;
            // Apply inter-particle damping proportional to radial velocity
            const float radialVel= (Vel[k0] - Vel[k1]).dot(distVecUnit);
            ForceMag[k0]+= (1.0f - distVal / maxForceLawRange) * D.UI[DampingRadRel___].F() * ForceLaw[Mat[k0]][0] * std::abs(radialVel) * surfArea;
            For[k0]-= (1.0f - distVal / maxForceLawRange) * D.UI[DampingRadRel___].F() * ForceLaw[Mat[k0]][0] * radialVel * surfArea * distVecUnit;
            // Increment neighbor count for display
            Neighbors[k0]++;
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
  const float particleMass= D.UI[MaterialDensity_].F() * std::pow(D.UI[LatticePitch____].F(), 3.0f);
  const float damping= std::max(D.UI[DampingVelRel___].F(), 0.0f);
  const bool useForceBC= D.UI[UseForceControl_].B();
  const int use2D= D.UI[ConstrainDim2D__].I();
  const Vec::Vec3<float> ForNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> ForPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());
  const Vec::Vec3<float> VelNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> VelPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());

  // Compute forces
  ComputeForces();                                        // f(t) = ∑ F(d)
  if (useForceBC)                                         // Check BC mode
    ApplyBCForces();                                      // Use force controller
  else                                                    // Check BC mode
    for (int k= 0; k < (int)Pos.size(); k++)              // Loop through elements
      if (BCFor[k] != 0)                                  // Check BC
        For[k]+= (BCFor[k] < 0) ? (ForNega) : (ForPosi);  // f(t)= f(t) + fext

  // Explicit forward Euler integration
  for (int k= 0; k < (int)Pos.size(); k++) {
    // Update velocities
    Vel[k]+= dt * For[k] / particleMass;           // v(t+1) = v(t) + Δt * f(t) / m
    if (damping > 0.0f)                            // Check damping mode
      Vel[k]= (1.0f - damping) * Vel[k];           // v(t+1)= ζ v(t+1)
    if (!useForceBC && BCVel[k] != 0)              // Check BC
      Vel[k]= (BCVel[k] < 0) ? VelNega : VelPosi;  // Overwrite velocity
    // Update positions
    Pos[k]+= dt * Vel[k];                          // x(t+1) = x(t) + Δt * v(t+1)
    if (!useForceBC && BCPos[k] != 0) {            // Check BC
      Pos[k]= Ref[k];                              // Overwrite position
      Vel[k]= Vec::Vec3<float>{0.0f, 0.0f, 0.0f};  // Reset velocity
    }  //
    else if (use2D > 0 && use2D <= 3) {      // Check 2D mode
      Pos[k][use2D - 1]= Ref[k][use2D - 1];  // Constrain to 2D
      Vel[k][use2D - 1]= 0.0f;               // Reset velocity
    }  //
  }

  // Advance time
  SimTime+= dt;

  // Scatter plot of sensor data
  if ((int)D.Scatter.size() < RunID) D.Scatter.resize(RunID);
  D.Scatter[RunID - 1].name= "ForceDisp";
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
  D.Scatter[RunID - 1].val.push_back(std::array<double, 2>{sumPos / (float)count, sumFor / (float)count});

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
  Field::Field3<double> tmpField(tmpNX, tmpNY, tmpNZ, 0.0);
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
            // TODO tweak decay function for cleaner results
            tmpField.at(x, y, z)+= std::max(0.0, 1.0 - (pos - Pos[k]).norm() / metaballSize);
          }
        }
      }
    }
  }

  MarchingCubes::ComputeMarchingCubes(tmpField.nX, tmpField.nY, tmpField.nZ, D.UI[MetaballIsoval__].F(), D.boxMin, D.boxMax, tmpField.data, Verts, Tris);

  if (D.UI[VerboseLevel____].I() >= 1) printf("MetaballT %f\n", Timer::PopTimer());
}
