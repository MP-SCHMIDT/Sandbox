#include "CompuFluidDyna.hpp"


// Standard lib
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
CompuFluidDyna::CompuFluidDyna() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void CompuFluidDyna::SetActiveProject() {
  if (!isActivProj) {
    D.UI.clear();
    D.UI.push_back(ParamUI("Scenario________", 0));      // Scenario ID, 0= load file, 1> hard coded scenarios
    D.UI.push_back(ParamUI("InputFile_______", 3));      // BMP file to load
    D.UI.push_back(ParamUI("ResolutionX_____", 1));      // Eulerian mesh resolution
    D.UI.push_back(ParamUI("ResolutionY_____", 100));    // Eulerian mesh resolution
    D.UI.push_back(ParamUI("ResolutionZ_____", 100));    // Eulerian mesh resolution
    D.UI.push_back(ParamUI("VoxelSize_______", 1.e-2));  // Element size
    D.UI.push_back(ParamUI("TimeStep________", 0.02));   // Simulation time step
    D.UI.push_back(ParamUI("SolvMaxIter_____", 32));     // Max number of solver iterations
    D.UI.push_back(ParamUI("SolvType________", 2));      // Flag to use Gauss Seidel (=0), Gradient Descent (=1) or Conjugate Gradient (=2)
    D.UI.push_back(ParamUI("SolvSOR_________", 1.8));    // Overrelaxation coefficient in Gauss Seidel solver
    D.UI.push_back(ParamUI("SolvTolRhs______", 0.0));    // Solver tolerance relative to RHS norm
    D.UI.push_back(ParamUI("SolvTolRel______", 1.e-3));  // Solver tolerance relative to initial guess
    D.UI.push_back(ParamUI("SolvTolAbs______", 0.0));    // Solver tolerance relative absolute value of residual magnitude
    D.UI.push_back(ParamUI("CoeffGravi______", 0.0));    // Magnitude of gravity in Z- direction
    D.UI.push_back(ParamUI("CoeffAdvec______", 5.0));    // 0= no advection, 1= linear advection, >1 MacCormack correction iterations
    D.UI.push_back(ParamUI("CoeffDiffuS_____", 1.e-4));  // Diffusion of smoke field, i.e. smoke spread/smear
    D.UI.push_back(ParamUI("CoeffDiffuV_____", 1.e-3));  // Diffusion of velocity field, i.e. viscosity
    D.UI.push_back(ParamUI("CoeffVorti______", 0.0));    // Vorticity confinement to avoid dissipation of energy in small scale vortices
    D.UI.push_back(ParamUI("CoeffProj_______", 1.0));    // Enable incompressibility projection
    D.UI.push_back(ParamUI("BCVelX__________", 0.0));    // Velocity value for voxels with enforced velocity
    D.UI.push_back(ParamUI("BCVelY__________", 1.0));    // Velocity value for voxels with enforced velocity
    D.UI.push_back(ParamUI("BCVelZ__________", 0.0));    // Velocity value for voxels with enforced velocity
    D.UI.push_back(ParamUI("BCPres__________", 1.0));    // Pressure value for voxels with enforced pressure
    D.UI.push_back(ParamUI("BCSmok__________", 1.0));    // Smoke value for voxels with enforced smoke
    D.UI.push_back(ParamUI("BCSmokTime______", 1.0));    // Period duration for input smoke oscillation
    D.UI.push_back(ParamUI("ObjectPosX______", 0.5));    // Coordinates for objects in hard coded scenarios
    D.UI.push_back(ParamUI("ObjectPosY______", 0.25));   // Coordinates for objects in hard coded scenarios
    D.UI.push_back(ParamUI("ObjectPosZ______", 0.5));    // Coordinates for objects in hard coded scenarios
    D.UI.push_back(ParamUI("ObjectSize0_____", 0.08));   // Size for objects in hard coded scenarios
    D.UI.push_back(ParamUI("ObjectSize1_____", 0.08));   // Size for objects in hard coded scenarios
    D.UI.push_back(ParamUI("ScaleFactor_____", 1.0));    // Scale factor for drawn geometry
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));    // Color factor for drawn geometry
    D.UI.push_back(ParamUI("ColorThresh_____", 0.0));    // Color cutoff drawn geometry
    D.UI.push_back(ParamUI("ColorMode_______", 1));      // Selector for the scalar field to be drawn
    D.UI.push_back(ParamUI("SliceDim________", 0));      // Enable model slicing along a dimension
    D.UI.push_back(ParamUI("SlicePlotX______", 0.5));    // Positions for the slices
    D.UI.push_back(ParamUI("SlicePlotY______", 0.5));    // Positions for the slices
    D.UI.push_back(ParamUI("SlicePlotZ______", 0.5));    // Positions for the slices
    D.UI.push_back(ParamUI("VerboseSolv_____", -0.5));   // Verbose mode for linear solvers
    D.UI.push_back(ParamUI("VerboseTime_____", -0.5));   // Verbose mode for linear solvers
    D.UI.push_back(ParamUI("VerboseLevel____", 0.0));    // Verbose mode
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
bool CompuFluidDyna::CheckAlloc() {
  if (D.UI[Scenario________].hasChanged()) isAllocated= false;
  if (D.UI[InputFile_______].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionX_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionY_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionZ_____].hasChanged()) isAllocated= false;
  if (D.UI[VoxelSize_______].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool CompuFluidDyna::CheckRefresh() {
  if (D.UI[BCVelX__________].hasChanged()) isRefreshed= false;
  if (D.UI[BCVelY__________].hasChanged()) isRefreshed= false;
  if (D.UI[BCVelZ__________].hasChanged()) isRefreshed= false;
  if (D.UI[BCPres__________].hasChanged()) isRefreshed= false;
  if (D.UI[BCSmok__________].hasChanged()) isRefreshed= false;
  if (D.UI[ObjectPosX______].hasChanged()) isRefreshed= false;
  if (D.UI[ObjectPosY______].hasChanged()) isRefreshed= false;
  if (D.UI[ObjectPosZ______].hasChanged()) isRefreshed= false;
  if (D.UI[ObjectSize0_____].hasChanged()) isRefreshed= false;
  if (D.UI[ObjectSize1_____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void CompuFluidDyna::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  nX= std::max(D.UI[ResolutionX_____].I(), 1);
  nY= std::max(D.UI[ResolutionY_____].I(), 1);
  nZ= std::max(D.UI[ResolutionZ_____].I(), 1);
  voxSize= std::max(D.UI[VoxelSize_______].F(), 1.e-6f);
  D.boxMin= {0.5f - 0.5f * (float)nX * voxSize, 0.5f - 0.5f * (float)nY * voxSize, 0.5f - 0.5f * (float)nZ * voxSize};
  D.boxMax= {0.5f + 0.5f * (float)nX * voxSize, 0.5f + 0.5f * (float)nY * voxSize, 0.5f + 0.5f * (float)nZ * voxSize};

  fluidDensity= 1.0f;

  // Allocate data
  Solid= Field::AllocField3D(nX, nY, nZ, false);
  VelBC= Field::AllocField3D(nX, nY, nZ, false);
  PreBC= Field::AllocField3D(nX, nY, nZ, false);
  SmoBC= Field::AllocField3D(nX, nY, nZ, false);
  VelXForced= Field::AllocField3D(nX, nY, nZ, 0.0f);
  VelYForced= Field::AllocField3D(nX, nY, nZ, 0.0f);
  VelZForced= Field::AllocField3D(nX, nY, nZ, 0.0f);
  PresForced= Field::AllocField3D(nX, nY, nZ, 0.0f);
  SmokForced= Field::AllocField3D(nX, nY, nZ, 0.0f);

  Dum0= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Dum1= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Dum2= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Dum3= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Dum4= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Vort= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Pres= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Dive= Field::AllocField3D(nX, nY, nZ, 0.0f);
  Smok= Field::AllocField3D(nX, nY, nZ, 0.0f);
  VelX= Field::AllocField3D(nX, nY, nZ, 0.0f);
  VelY= Field::AllocField3D(nX, nY, nZ, 0.0f);
  VelZ= Field::AllocField3D(nX, nY, nZ, 0.0f);
  CurX= Field::AllocField3D(nX, nY, nZ, 0.0f);
  CurY= Field::AllocField3D(nX, nY, nZ, 0.0f);
  CurZ= Field::AllocField3D(nX, nY, nZ, 0.0f);
  AdvX= Field::AllocField3D(nX, nY, nZ, 0.0f);
  AdvY= Field::AllocField3D(nX, nY, nZ, 0.0f);
  AdvZ= Field::AllocField3D(nX, nY, nZ, 0.0f);
}


// Refresh the project
void CompuFluidDyna::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Initialize scenario values
  simTime= 0;
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        Solid[x][y][z]= VelBC[x][y][z]= PreBC[x][y][z]= SmoBC[x][y][z]= false;
        VelXForced[x][y][z]= VelYForced[x][y][z]= VelZForced[x][y][z]= PresForced[x][y][z]= SmokForced[x][y][z]= 0.0f;
      }
    }
  }

  CompuFluidDyna::InitializeScenario();

  // Enforce scenario validity
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z]) {
          VelBC[x][y][z]= PreBC[x][y][z]= SmoBC[x][y][z]= false;
          VelXForced[x][y][z]= VelYForced[x][y][z]= VelZForced[x][y][z]= PresForced[x][y][z]= SmokForced[x][y][z]= 0.0f;
        }
      }
    }
  }

  // Apply BC on fields
  ApplyBC(FieldID::IDSmok, Smok);
  ApplyBC(FieldID::IDVelX, VelX);
  ApplyBC(FieldID::IDVelY, VelY);
  ApplyBC(FieldID::IDVelZ, VelZ);
  ApplyBC(FieldID::IDPres, Dive);
  ApplyBC(FieldID::IDPres, Pres);
}


// Handle keypress
void CompuFluidDyna::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Animate the project
void CompuFluidDyna::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Get simulation parameters
  const int maxIter= std::max(D.UI[SolvMaxIter_____].I(), 0);
  const float timestep= D.UI[TimeStep________].F();
  const float coeffDiffu= std::max(D.UI[CoeffDiffuS_____].F(), 0.0f);
  const float coeffVisco= std::max(D.UI[CoeffDiffuV_____].F(), 0.0f);
  const float coeffVorti= D.UI[CoeffVorti______].F();
  simTime+= D.UI[TimeStep________].F();

  // Update periodic smoke in inlet
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  ApplyBC(FieldID::IDSmok, Smok);
  if (D.UI[VerboseTime_____].B()) printf("%f T ApplyBC\n", Timer::PopTimer());

  // Incompressible Navier Stokes
  // ∂vel/∂t = - (vel · ∇) vel + visco ∇²vel − 1/ρ ∇press + f
  // ∇ · vel = 0

  // Advection steps
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  if (D.UI[CoeffAdvec______].B()) {
    AdvectField(FieldID::IDSmok, timestep, VelX, VelY, VelZ, Smok);
  }
  if (D.UI[CoeffAdvec______].B()) {
    std::vector<std::vector<std::vector<float>>> oldVelX= VelX;
    std::vector<std::vector<std::vector<float>>> oldVelY= VelY;
    std::vector<std::vector<std::vector<float>>> oldVelZ= VelZ;
    if (nX > 1) AdvectField(FieldID::IDVelX, timestep, oldVelX, oldVelY, oldVelZ, VelX);
    if (nY > 1) AdvectField(FieldID::IDVelY, timestep, oldVelX, oldVelY, oldVelZ, VelY);
    if (nZ > 1) AdvectField(FieldID::IDVelZ, timestep, oldVelX, oldVelY, oldVelZ, VelZ);
  }
  if (D.UI[VerboseTime_____].B()) printf("%f T AdvectField\n", Timer::PopTimer());

  // Diffusion steps
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  if (D.UI[CoeffDiffuS_____].B()) {
    // (Id - diffu Δt ∇²) smo = smo
    std::vector<std::vector<std::vector<float>>> oldSmoke= Smok;
    if (D.UI[SolvType________].I() == 0) {
      GaussSeidelSolve(FieldID::IDSmok, maxIter, timestep, true, coeffDiffu, oldSmoke, Smok);
    }
    else if (D.UI[SolvType________].I() == 1) {
      GradientDescentSolve(FieldID::IDSmok, maxIter, timestep, true, coeffDiffu, oldSmoke, Smok);
    }
    else {
      ConjugateGradientSolve(FieldID::IDSmok, maxIter, timestep, true, coeffDiffu, oldSmoke, Smok);
    }
  }
  if (D.UI[CoeffDiffuV_____].B()) {
    // (Id - visco Δt ∇²) vel = vel
    std::vector<std::vector<std::vector<float>>> oldVelX= VelX;
    std::vector<std::vector<std::vector<float>>> oldVelY= VelY;
    std::vector<std::vector<std::vector<float>>> oldVelZ= VelZ;
    if (D.UI[SolvType________].I() == 0) {
      if (nX > 1) GaussSeidelSolve(FieldID::IDVelX, maxIter, timestep, true, coeffVisco, oldVelX, VelX);
      if (nY > 1) GaussSeidelSolve(FieldID::IDVelY, maxIter, timestep, true, coeffVisco, oldVelY, VelY);
      if (nZ > 1) GaussSeidelSolve(FieldID::IDVelZ, maxIter, timestep, true, coeffVisco, oldVelZ, VelZ);
    }
    else if (D.UI[SolvType________].I() == 1) {
      if (nX > 1) GradientDescentSolve(FieldID::IDVelX, maxIter, timestep, true, coeffVisco, oldVelX, VelX);
      if (nY > 1) GradientDescentSolve(FieldID::IDVelY, maxIter, timestep, true, coeffVisco, oldVelY, VelY);
      if (nZ > 1) GradientDescentSolve(FieldID::IDVelZ, maxIter, timestep, true, coeffVisco, oldVelZ, VelZ);
    }
    else {
      if (nX > 1) ConjugateGradientSolve(FieldID::IDVelX, maxIter, timestep, true, coeffVisco, oldVelX, VelX);
      if (nY > 1) ConjugateGradientSolve(FieldID::IDVelY, maxIter, timestep, true, coeffVisco, oldVelY, VelY);
      if (nZ > 1) ConjugateGradientSolve(FieldID::IDVelZ, maxIter, timestep, true, coeffVisco, oldVelZ, VelZ);
    }
  }
  if (D.UI[VerboseTime_____].B()) printf("%f T Diffusion\n", Timer::PopTimer());

  // Vorticity step
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  if (D.UI[CoeffVorti______].B()) {
    VorticityConfinement(timestep, coeffVorti, VelX, VelY, VelZ);
  }
  if (D.UI[VerboseTime_____].B()) printf("%f T VorticityConfinement\n", Timer::PopTimer());

  // External forces
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  if (D.UI[CoeffGravi______].B()) {
    ExternalForces();
  }
  if (D.UI[VerboseTime_____].B()) printf("%f T ExternalForces\n", Timer::PopTimer());

  // Projection step
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  if (D.UI[CoeffProj_______].B()) {
    ProjectField(maxIter, timestep, VelX, VelY, VelZ);
  }
  if (D.UI[VerboseTime_____].B()) printf("%f T ProjectField\n", Timer::PopTimer());

  // Compute field data for display
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  ComputeVelocityDivergence();
  if (D.UI[VerboseTime_____].B()) printf("%f T Divergence\n", Timer::PopTimer());
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  ComputeVelocityCurlVorticity();
  if (D.UI[VerboseTime_____].B()) printf("%f T CurlVorticity\n", Timer::PopTimer());

  // TODO Compute fluid density to check if constant as it should be in incompressible case

  // TODO Test heuristic optimization of solid regions
  // https://open-research-europe.ec.europa.eu/articles/3-156

  // TODO Introduce solid interface normals calculations to better handle BC on sloped geometry ?
  // TODO Test with flow separation scenarios ?

  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  CompuFluidDyna::SetUpUIData();
  if (D.UI[VerboseTime_____].B()) printf("%f T SetUpUIData\n", Timer::PopTimer());

  if (D.UI[VerboseTime_____].B()) printf("\n");
}


// Draw the project
void CompuFluidDyna::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the voxels
  if (D.displayMode1) {
    glLineWidth(2.0f);
    // Set the scene transformation
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    // Sweep the field
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (D.UI[SliceDim________].I() == 1 && x != (int)std::round(D.UI[SlicePlotX______].F() * nX)) continue;
          if (D.UI[SliceDim________].I() == 2 && y != (int)std::round(D.UI[SlicePlotY______].F() * nY)) continue;
          if (D.UI[SliceDim________].I() == 3 && z != (int)std::round(D.UI[SlicePlotZ______].F() * nZ)) continue;
          // Set the voxel color components
          float r= 0.4f, g= 0.4f, b= 0.4f;
          if (PreBC[x][y][z]) r= 0.7f;
          if (VelBC[x][y][z]) g= 0.7f;
          if (SmoBC[x][y][z]) b= 0.7f;
          // Draw the cube
          if (Solid[x][y][z] || PreBC[x][y][z] || VelBC[x][y][z] || SmoBC[x][y][z]) {
            glColor3f(r, g, b);
            glPushMatrix();
            glTranslatef((float)x, (float)y, (float)z);
            glutWireCube(1.0);
            glPopMatrix();
          }
        }
      }
    }
    glPopMatrix();
    glLineWidth(1.0f);
  }

  // Draw the scalar fields
  if (D.displayMode2) {
    // Set the scene transformation
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    if (nX == 1) glScalef(0.1f, 1.0f, 1.0f);
    if (nY == 1) glScalef(1.0f, 0.1f, 1.0f);
    if (nZ == 1) glScalef(1.0f, 1.0f, 0.1f);
    // Sweep the field
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (D.UI[SliceDim________].I() == 1 && x != (int)std::round(D.UI[SlicePlotX______].F() * nX)) continue;
          if (D.UI[SliceDim________].I() == 2 && y != (int)std::round(D.UI[SlicePlotY______].F() * nY)) continue;
          if (D.UI[SliceDim________].I() == 3 && z != (int)std::round(D.UI[SlicePlotZ______].F() * nZ)) continue;
          if (Solid[x][y][z] && D.UI[ColorThresh_____].F() == 0.0) continue;
          float r= 0.0f, g= 0.0f, b= 0.0f;
          // Color by smoke
          if (D.UI[ColorMode_______].I() == 1) {
            if (std::abs(Smok[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToPlasma(0.5f + 0.5f * Smok[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by velocity magnitude
          if (D.UI[ColorMode_______].I() == 2) {
            Vec::Vec3<float> vec(VelX[x][y][z], VelY[x][y][z], VelZ[x][y][z]);
            if (vec.norm() < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToJetBrightSmooth(vec.norm() * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by pressure
          if (D.UI[ColorMode_______].I() == 3) {
            if (std::abs(Pres[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToBlueToRed(0.5f + 0.5f * Pres[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by divergence
          if (D.UI[ColorMode_______].I() == 4) {
            if (std::abs(Dive[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToGreenToRed(0.5f + 0.5f * Dive[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by vorticity
          if (D.UI[ColorMode_______].I() == 5) {
            if (std::abs(Vort[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToJetBrightSmooth(Vort[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by vel in X
          if (D.UI[ColorMode_______].I() == 6) {
            if (std::abs(VelX[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToJetBrightSmooth(0.5f + 0.5f * VelX[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by vel in Y
          if (D.UI[ColorMode_______].I() == 7) {
            if (std::abs(VelY[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToJetBrightSmooth(0.5f + 0.5f * VelY[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by vel in Z
          if (D.UI[ColorMode_______].I() == 8) {
            if (std::abs(VelZ[x][y][z]) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToJetBrightSmooth(0.5f + 0.5f * VelZ[x][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          }
          // Color by dummy values
          if (D.UI[ColorMode_______].I() >= 10 && D.UI[ColorMode_______].I() <= 14) {
            float val= 0.0f;
            if (D.UI[ColorMode_______].I() == 10) val= Dum0[x][y][z];
            if (D.UI[ColorMode_______].I() == 11) val= Dum1[x][y][z];
            if (D.UI[ColorMode_______].I() == 12) val= Dum2[x][y][z];
            if (D.UI[ColorMode_______].I() == 13) val= Dum3[x][y][z];
            if (D.UI[ColorMode_______].I() == 14) val= Dum4[x][y][z];
            if (std::abs(val) < D.UI[ColorThresh_____].F()) continue;
            Colormap::RatioToJetBrightSmooth(0.5f + 0.5f * val * D.UI[ColorFactor_____].F(), r, g, b);
          }
          glColor3f(r, g, b);
          glPushMatrix();
          glTranslatef((float)x, (float)y, (float)z);
          glutSolidCube(1.0);
          glPopMatrix();
        }
      }
    }
    glPopMatrix();
  }

  // Draw the velocity field
  if (D.displayMode3) {
    // Set the scene transformation
    glPushMatrix();
    if (nX == 1) glTranslatef(voxSize, 0.0f, 0.0f);
    if (nY == 1) glTranslatef(0.0f, voxSize, 0.0f);
    if (nZ == 1) glTranslatef(0.0f, 0.0f, voxSize);
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    // Sweep the field
    constexpr int nbLineWidths= 3;
    for (int k= 0; k < nbLineWidths; k++) {
      const float segmentRelLength= 1.0f - (float)k / (float)nbLineWidths;
      glLineWidth((float)k + 1.0f);
      glBegin(GL_LINES);
      for (int x= 0; x < nX; x++) {
        for (int y= 0; y < nY; y++) {
          for (int z= 0; z < nZ; z++) {
            if (D.UI[SliceDim________].I() == 1 && x != (int)std::round(D.UI[SlicePlotX______].F() * nX)) continue;
            if (D.UI[SliceDim________].I() == 2 && y != (int)std::round(D.UI[SlicePlotY______].F() * nY)) continue;
            if (D.UI[SliceDim________].I() == 3 && z != (int)std::round(D.UI[SlicePlotZ______].F() * nZ)) continue;
            if (Solid[x][y][z] && D.UI[ColorThresh_____].F() == 0.0) continue;
            // Draw the velocity field
            Vec::Vec3<float> vec(VelX[x][y][z], VelY[x][y][z], VelZ[x][y][z]);
            if (std::abs(D.UI[SliceDim________].I()) == 1) vec[0]= 0.0f;
            if (std::abs(D.UI[SliceDim________].I()) == 2) vec[1]= 0.0f;
            if (std::abs(D.UI[SliceDim________].I()) == 3) vec[2]= 0.0f;
            if (vec.normSquared() > 0.0f) {
              float r= 0.0f, g= 0.0f, b= 0.0f;
              Colormap::RatioToJetBrightSmooth(vec.norm() * D.UI[ColorFactor_____].F(), r, g, b);
              glColor3f(r, g, b);
              Vec::Vec3<float> pos((float)x, (float)y, (float)z);
              glVertex3fv(pos.array());
              // glVertex3fv(pos + vec * segmentRelLength * D.UI[ScaleFactor_____].F());
              // glVertex3fv(pos + vec.normalized() * segmentRelLength * D.UI[ScaleFactor_____].F());
              // glVertex3fv(pos + vec.normalized() * segmentRelLength * D.UI[ScaleFactor_____].F() * std::log(vec.norm() + 1.0f));
              glVertex3fv(pos + vec.normalized() * segmentRelLength * D.UI[ScaleFactor_____].F() * std::sqrt(vec.norm()));
            }
          }
        }
      }
      glEnd();
    }
    glLineWidth(1.0f);
    glPopMatrix();
  }

  // Draw the advection source field
  if (D.displayMode4) {
    // Set the scene transformation
    glPushMatrix();
    if (nX == 1) glTranslatef(voxSize, 0.0f, 0.0f);
    if (nY == 1) glTranslatef(0.0f, voxSize, 0.0f);
    if (nZ == 1) glTranslatef(0.0f, 0.0f, voxSize);
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    // Sweep the field
    glBegin(GL_LINES);
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (D.UI[SliceDim________].I() == 1 && x != (int)std::round(D.UI[SlicePlotX______].F() * nX)) continue;
          if (D.UI[SliceDim________].I() == 2 && y != (int)std::round(D.UI[SlicePlotY______].F() * nY)) continue;
          if (D.UI[SliceDim________].I() == 3 && z != (int)std::round(D.UI[SlicePlotZ______].F() * nZ)) continue;
          // Draw the velocity field
          Vec::Vec3<float> vec(AdvX[x][y][z], AdvY[x][y][z], AdvZ[x][y][z]);
          if (std::abs(D.UI[SliceDim________].I()) == 1) vec[0]= 0.0f;
          if (std::abs(D.UI[SliceDim________].I()) == 2) vec[1]= 0.0f;
          if (std::abs(D.UI[SliceDim________].I()) == 3) vec[2]= 0.0f;
          if (vec.normSquared() > 0.0f) {
            const float r= 0.5f - vec[0];
            const float g= 0.5f - vec[1];
            const float b= 0.5f - vec[2];
            glColor3f(r, g, b);
            Vec::Vec3<float> pos((float)x, (float)y, (float)z);
            glVertex3fv(pos.array());
            glVertex3fv(pos + vec);
          }
        }
      }
    }
    glEnd();
    glPopMatrix();
  }
}
