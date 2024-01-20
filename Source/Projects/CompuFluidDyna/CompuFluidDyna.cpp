#include "CompuFluidDyna.hpp"


// Standard lib
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
    D.UI.push_back(ParamUI("CoeffAdvec______", 5));      // 0= no advection, 1= linear advection, >1 MacCormack correction iterations
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
    D.UI.push_back(ParamUI("VerboseLevel____", 0));      // Verbose mode
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

  // Step forward in the fluid simulation
  RunSimulationStep();

  // Update the plots
  if (D.UI[VerboseTime_____].B()) Timer::PushTimer();
  CompuFluidDyna::UpdateUIData();
  if (D.UI[VerboseTime_____].B()) printf("%f T UpdateUIData\n", Timer::PopTimer());

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

void CompuFluidDyna::UpdateUIData() {
  // Draw the scatter data
  const int yCursor= std::min(std::max((int)std::round((float)(nY - 1) * D.UI[SlicePlotY______].F()), 0), nY - 1);
  const int zCursor= std::min(std::max((int)std::round((float)(nZ - 1) * D.UI[SlicePlotZ______].F()), 0), nZ - 1);
  D.scatLegend.resize(4);
  D.scatLegend[0]= "Horiz VZ";
  D.scatLegend[1]= "Verti VY";
  D.scatLegend[2]= "Horiz P";
  D.scatLegend[3]= "Verti P";
  D.scatData.resize(4);
  for (int k= 0; k < (int)D.scatData.size(); k++)
    D.scatData[k].clear();
  if (nZ > 1) {
    for (int y= 0; y < nY; y++) {
      D.scatData[0].push_back(std::array<double, 2>({(double)y / (double)(nY - 1), VelZ[nX / 2][y][zCursor]}));
      D.scatData[2].push_back(std::array<double, 2>({(double)y / (double)(nY - 1), Pres[nX / 2][y][zCursor]}));
    }
  }
  if (nY > 1) {
    for (int z= 0; z < nZ; z++) {
      D.scatData[1].push_back(std::array<double, 2>({VelY[nX / 2][yCursor][z], (double)z / (double)(nZ - 1)}));
      D.scatData[3].push_back(std::array<double, 2>({Pres[nX / 2][yCursor][z], (double)z / (double)(nZ - 1)}));
    }
  }

  // Add hard coded experimental values for lid driven cavity flow benchmark
  if (D.UI[Scenario________].I() == 3) {
    // Clear unnecessary scatter data
    D.scatData[2].clear();
    D.scatData[3].clear();
    // Allocate required scatter data
    D.scatLegend.resize(8);
    D.scatLegend[4]= "Ghia Re1k";
    D.scatLegend[5]= "Ghia Re1k";
    D.scatLegend[6]= "Ertu Re1k";
    D.scatLegend[7]= "Ertu Re1k";
    D.scatData.resize(8);
    D.scatData[4].clear();
    D.scatData[5].clear();
    D.scatData[6].clear();
    D.scatData[7].clear();
    // TODO add vorticity data from https://www.acenumerics.com/the-benchmarks.html
    // Data from Ghia 1982 http://www.msaidi.ir/upload/Ghia1982.pdf
    const std::vector<double> GhiaData0X({0.0000, +0.0625, +0.0703, +0.0781, +0.0938, +0.1563, +0.2266, +0.2344, +0.5000, +0.8047, +0.8594, +0.9063, +0.9453, +0.9531, +0.9609, +0.9688, +1.0000});  // coord on horiz slice
    const std::vector<double> GhiaData1Y({0.0000, +0.0547, +0.0625, +0.0703, +0.1016, +0.1719, +0.2813, +0.4531, +0.5000, +0.6172, +0.7344, +0.8516, +0.9531, +0.9609, +0.9688, +0.9766, +1.0000});  // coord on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.0923, +0.1009, +0.1089, +0.1232, +0.1608, +0.1751, +0.1753, +0.0545, -0.2453, -0.2245, -0.1691, -0.1031, -0.0886, -0.0739, -0.0591, +0.0000});  // Re 100   verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.0372, -0.0419, -0.0478, -0.0643, -0.1015, -0.1566, -0.2109, -0.2058, -0.1364, +0.0033, +0.2315, +0.6872, +0.7372, +0.7887, +0.8412, +1.0000});  // Re 100   horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.1836, +0.1971, +0.2092, +0.2297, +0.2812, +0.3020, +0.3017, +0.0519, -0.3860, -0.4499, -0.2383, -0.2285, -0.1925, -0.1566, -0.1215, +0.0000});  // Re 400   verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.0819, -0.0927, -0.1034, -0.1461, -0.2430, -0.3273, -0.1712, -0.1148, +0.0214, +0.1626, +0.2909, +0.5589, +0.6176, +0.6844, +0.7584, +1.0000});  // Re 400   horiz vel on verti slice
    const std::vector<double> GhiaData0Y({0.0000, +0.2749, +0.2901, +0.3035, +0.3263, +0.3710, +0.3308, +0.3224, +0.0253, -0.3197, -0.4267, -0.5150, -0.3919, -0.3371, -0.2767, -0.2139, +0.0000});  // Re 1000  verti vel on horiz slice
    const std::vector<double> GhiaData1X({0.0000, -0.1811, -0.2020, -0.2222, -0.2973, -0.3829, -0.2781, -0.1065, -0.0608, +0.0570, +0.1872, +0.3330, +0.4660, +0.5112, +0.5749, +0.6593, +1.0000});  // Re 1000  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.3956, +0.4092, +0.4191, +0.4277, +0.3712, +0.2903, +0.2819, +0.0100, -0.3118, -0.3740, -0.4431, -0.5405, -0.5236, -0.4743, -0.3902, +0.0000});  // Re 3200  verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.3241, -0.3534, -0.3783, -0.4193, -0.3432, -0.2443, -0.8664, -0.0427, +0.0716, +0.1979, +0.3468, +0.4610, +0.4655, +0.4830, +0.5324, +1.0000});  // Re 3200  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.4245, +0.4333, +0.4365, +0.4295, +0.3537, +0.2807, +0.2728, +0.0095, -0.3002, -0.3621, -0.4144, -0.5288, -0.5541, -0.5507, -0.4977, +0.0000});  // Re 5000  verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.4117, -0.4290, -0.4364, -0.4044, -0.3305, -0.2286, -0.0740, -0.0304, +0.0818, +0.2009, +0.3356, +0.4604, +0.4599, +0.4612, +0.4822, +1.0000});  // Re 5000  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.4398, +0.4403, +0.4356, +0.4182, +0.3506, +0.2812, +0.2735, +0.0082, -0.3045, -0.3621, -0.4105, -0.4859, -0.5235, -0.5522, -0.5386, +0.0000});  // Re 7500  verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.4315, -0.4359, -0.4303, -0.3832, -0.3239, -0.2318, -0.0750, -0.0380, +0.0834, +0.2059, +0.3423, +0.4717, +0.4732, +0.4705, +0.4724, +1.0000});  // Re 7500  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.4398, +0.4373, +0.4312, +0.4149, +0.3507, +0.2800, +0.2722, +0.0083, -0.3072, -0.3674, -0.4150, -0.4586, -0.4910, -0.5299, -0.5430, +0.0000});  // Re 10000 verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.4274, -0.4254, -0.4166, -0.3800, -0.3271, -0.2319, -0.0754, +0.0311, +0.0834, +0.2067, +0.3464, +0.4780, +0.4807, +0.4778, +0.4722, +1.0000});  // Re 10000 horiz vel on verti slice
    // Data from Erturk 2005 https://arxiv.org/pdf/physics/0505121.pdf
    const std::vector<double> ErtuData0X({0.0000, +0.0150, +0.0300, +0.0450, +0.0600, +0.0750, +0.0900, +0.1050, +0.1200, +0.1350, +0.1500, +0.5000, +0.8500, +0.8650, +0.8800, +0.8950, +0.9100, +0.9250, +0.9400, +0.9550, +0.9700, +0.9850, +1.0000});  // coord on horiz slice
    const std::vector<double> ErtuData1Y({0.0000, +0.0200, +0.0400, +0.0600, +0.0800, +0.1000, +0.1200, +0.1400, +0.1600, +0.1800, +0.2000, +0.5000, +0.9000, +0.9100, +0.9200, +0.9300, +0.9400, +0.9500, +0.9600, +0.9700, +0.9800, +0.9900, +1.0000});  // coord on verti slice
    const std::vector<double> ErtuData0Y({0.0000, +0.1019, +0.1792, +0.2349, +0.2746, +0.3041, +0.3273, +0.3460, +0.3605, +0.3705, +0.3756, +0.0258, -0.4028, -0.4407, -0.4803, -0.5132, -0.5263, -0.5052, -0.4417, -0.3400, -0.2173, -0.0973, +0.0000});  // Re 1000  verti vel on horiz slice
    const std::vector<double> ErtuData1X({0.0000, -0.0757, -0.1392, -0.1951, -0.2472, -0.2960, -0.3381, -0.3690, -0.3854, -0.3869, -0.3756, -0.0620, +0.3838, +0.3913, +0.3993, +0.4101, +0.4276, +0.4582, +0.5102, +0.5917, +0.7065, +0.8486, +1.0000});  // Re 1000  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.1607, +0.2633, +0.3238, +0.3649, +0.3950, +0.4142, +0.4217, +0.4187, +0.4078, +0.3918, +0.0160, -0.3671, -0.3843, -0.4042, -0.4321, -0.4741, -0.5268, -0.5603, -0.5192, -0.3725, -0.1675, +0.0000});  // Re 2500  verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.1517, -0.2547, -0.3372, -0.3979, -0.4250, -0.4200, -0.3965, -0.3688, -0.3439, -0.3228, -0.0403, +0.4141, +0.4256, +0.4353, +0.4424, +0.4470, +0.4506, +0.4607, +0.4971, +0.5924, +0.7704, +1.0000});  // Re 2500  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2160, +0.3263, +0.3868, +0.4258, +0.4426, +0.4403, +0.4260, +0.4070, +0.3878, +0.3699, +0.0117, -0.3624, -0.3806, -0.3982, -0.4147, -0.4318, -0.4595, -0.5139, -0.5700, -0.5019, -0.2441, +0.0000});  // Re 5000  verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.2223, -0.3480, -0.4272, -0.4419, -0.4168, -0.3876, -0.3652, -0.3467, -0.3285, -0.3100, -0.0319, +0.4155, +0.4307, +0.4452, +0.4582, +0.4683, +0.4738, +0.4739, +0.4749, +0.5159, +0.6866, +1.0000});  // Re 5000  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2509, +0.3608, +0.4210, +0.4494, +0.4495, +0.4337, +0.4137, +0.3950, +0.3779, +0.3616, +0.0099, -0.3574, -0.3755, -0.3938, -0.4118, -0.4283, -0.4443, -0.4748, -0.5434, -0.5550, -0.2991, +0.0000});  // Re 7500  verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.2633, -0.3980, -0.4491, -0.4284, -0.3978, -0.3766, -0.3587, -0.3406, -0.3222, -0.3038, -0.0287, +0.4123, +0.4275, +0.4431, +0.4585, +0.4723, +0.4824, +0.4860, +0.4817, +0.4907, +0.6300, +1.0000});  // Re 7500  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2756, +0.3844, +0.4409, +0.4566, +0.4449, +0.4247, +0.4056, +0.3885, +0.3722, +0.3562, +0.0088, -0.3538, -0.3715, -0.3895, -0.4078, -0.4256, -0.4411, -0.4592, -0.5124, -0.5712, -0.3419, +0.0000});  // Re 10000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.2907, -0.4259, -0.4469, -0.4142, -0.3899, -0.3721, -0.3543, -0.3361, -0.3179, -0.2998, -0.0268, +0.4095, +0.4243, +0.4398, +0.4556, +0.4711, +0.4843, +0.4917, +0.4891, +0.4837, +0.5891, +1.0000});  // Re 10000 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2940, +0.4018, +0.4522, +0.4563, +0.4383, +0.4180, +0.4004, +0.3840, +0.3678, +0.3519, +0.0080, -0.3508, -0.3682, -0.3859, -0.4040, -0.4221, -0.4388, -0.4534, -0.4899, -0.5694, -0.3762, +0.0000});  // Re 12500 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3113, -0.4407, -0.4380, -0.4054, -0.3859, -0.3685, -0.3506, -0.3326, -0.3146, -0.2967, -0.0256, +0.4070, +0.4216, +0.4366, +0.4523, +0.4684, +0.4833, +0.4937, +0.4941, +0.4833, +0.5587, +1.0000});  // Re 12500 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3083, +0.4152, +0.4580, +0.4529, +0.4323, +0.4132, +0.3964, +0.3801, +0.3641, +0.3483, +0.0074, -0.3481, -0.3654, -0.3828, -0.4005, -0.4186, -0.4361, -0.4505, -0.4754, -0.5593, -0.4041, +0.0000});  // Re 15000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3278, -0.4474, -0.4286, -0.4001, -0.3827, -0.3652, -0.3474, -0.3297, -0.3119, -0.2942, -0.0247, +0.4047, +0.4190, +0.4338, +0.4492, +0.4653, +0.4811, +0.4937, +0.4969, +0.4850, +0.5358, +1.0000});  // Re 15000 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3197, +0.4254, +0.4602, +0.4484, +0.4273, +0.4093, +0.3929, +0.3767, +0.3608, +0.3452, +0.0069, -0.3457, -0.3627, -0.3800, -0.3975, -0.4153, -0.4331, -0.4482, -0.4664, -0.5460, -0.4269, +0.0000});  // Re 17500 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3412, -0.4490, -0.4206, -0.3965, -0.3797, -0.3622, -0.3446, -0.3271, -0.3096, -0.2920, -0.0240, +0.4024, +0.4166, +0.4312, +0.4463, +0.4622, +0.4784, +0.4925, +0.4982, +0.4871, +0.5183, +1.0000});  // Re 17500 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3290, +0.4332, +0.4601, +0.4438, +0.4232, +0.4060, +0.3897, +0.3736, +0.3579, +0.3423, +0.0065, -0.3434, -0.3603, -0.3774, -0.3946, -0.4122, -0.4300, -0.4459, -0.4605, -0.5321, -0.4457, +0.0000});  // Re 20000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3523, -0.4475, -0.4143, -0.3936, -0.3769, -0.3595, -0.3422, -0.3248, -0.3074, -0.2899, -0.0234, +0.4001, +0.4142, +0.4287, +0.4436, +0.4592, +0.4754, +0.4906, +0.4985, +0.4889, +0.5048, +1.0000});  // Re 20000 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3323, +0.4357, +0.4596, +0.4420, +0.4218, +0.4048, +0.3885, +0.3725, +0.3567, +0.3413, +0.0063, -0.3425, -0.3593, -0.3764, -0.3936, -0.4110, -0.4287, -0.4449, -0.4588, -0.5266, -0.4522, +0.0000});  // Re 21000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3562, -0.4463, -0.4121, -0.3925, -0.3758, -0.3585, -0.3412, -0.3239, -0.3066, -0.2892, -0.0232, +0.3992, +0.4132, +0.4277, +0.4425, +0.4580, +0.4742, +0.4897, +0.4983, +0.4895, +0.5003, +1.0000});  // Re 21000 horiz vel on verti slice
    // Add hard coded experimental values in the scatter plot
    for (int k= 0; k < (int)GhiaData0X.size(); k++) {
      D.scatData[4].push_back(std::array<double, 2>({GhiaData0X[k], GhiaData0Y[k]}));
      D.scatData[5].push_back(std::array<double, 2>({GhiaData1X[k], GhiaData1Y[k]}));
    }
    for (int k= 0; k < (int)ErtuData0X.size(); k++) {
      D.scatData[6].push_back(std::array<double, 2>({ErtuData0X[k], ErtuData0Y[k]}));
      D.scatData[7].push_back(std::array<double, 2>({ErtuData1X[k], ErtuData1Y[k]}));
    }
  }

  // Add hard coded analytical values for Poiseuille flow benchmark
  if (D.UI[Scenario________].I() == 6) {
    // Clear unnecessary scatter data
    D.scatData[0].clear();
    D.scatData[3].clear();
    // Allocate required scatter data
    D.scatLegend.resize(6);
    D.scatLegend[4]= "Analy VY";
    D.scatLegend[5]= "Analy P";
    D.scatData.resize(6);
    D.scatData[4].clear();
    D.scatData[5].clear();
    // Add analytical values in the scatter plot
    const float press0= D.UI[BCPres__________].F();
    const float press1= -D.UI[BCPres__________].F();
    const float kinVisco= D.UI[CoeffDiffuV_____].F();
    if (nY > 1) {
      for (int z= 0; z < nZ; z++) {
        const float width= voxSize * (float)(nY - 1);
        const float height= voxSize * (float)(nZ - 1);
        const float posZ= (float)z * voxSize;
        const float pressDiff= (press1 - press0) / width;
        const float analyVelY= -pressDiff * (1.0f / (2.0f * kinVisco)) * posZ * (height - posZ);
        D.scatData[4].push_back(std::array<double, 2>({analyVelY, (double)z / (double)(nZ - 1)}));
      }
    }
    if (nZ > 1) {
      for (int y= 0; y < nY; y++) {
        const float analyP= press0 + (press1 - press0) * (float)y / (float)(nY - 1);
        D.scatData[5].push_back(std::array<double, 2>({(double)y / (double)(nY - 1), analyP}));
      }
    }
  }

  if (!D.UI[VerboseSolv_____].B()) {
    // Draw the plot data
    D.plotData.resize(5);
    D.plotLegend.resize(5);
    D.plotLegend[0]= "VelMag";
    D.plotLegend[1]= "Smok";
    D.plotLegend[2]= "Pres";
    D.plotLegend[3]= "DiveAbs";
    D.plotLegend[4]= "Vorti";
    if (D.plotData[0].size() < 1000) {
      for (int k= 0; k < (int)D.plotLegend.size(); k++)
        D.plotData[k].push_back(0.0f);
      for (int x= 0; x < nX; x++) {
        for (int y= 0; y < nY; y++) {
          for (int z= 0; z < nZ; z++) {
            D.plotData[0][D.plotData[0].size() - 1]+= std::sqrt(VelX[x][y][z] * VelX[x][y][z] + VelY[x][y][z] * VelY[x][y][z] + VelZ[x][y][z] * VelZ[x][y][z]);
            D.plotData[1][D.plotData[1].size() - 1]+= Smok[x][y][z];
            D.plotData[2][D.plotData[2].size() - 1]+= Pres[x][y][z];
            D.plotData[3][D.plotData[3].size() - 1]+= std::abs(Dive[x][y][z]);
            D.plotData[4][D.plotData[4].size() - 1]+= Vort[x][y][z];
          }
        }
      }
    }
  }
}


void CompuFluidDyna::InitializeScenario() {
  // Get scenario ID and optionnally load bitmap file
  const int scenarioType= D.UI[Scenario________].I();
  const int inputFile= D.UI[InputFile_______].I();
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (scenarioType == 0) {
    if (inputFile == 0) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/NACA.bmp", imageRGBA, false);
    if (inputFile == 1) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/Nozzle.bmp", imageRGBA, false);
    if (inputFile == 2) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/Pipe.bmp", imageRGBA, false);
    if (inputFile == 3) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/PipePressureBC.bmp", imageRGBA, false);
    if (inputFile == 4) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/TeslaValve.bmp", imageRGBA, false);
  }

  // Set scenario values
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Scenario from loaded BMP file
        if (scenarioType == 0 && !imageRGBA.empty()) {
          const float posW= (float)(imageRGBA.size() - 1) * ((float)y + 0.5f) / (float)nY;
          const float posH= (float)(imageRGBA[0].size() - 1) * ((float)z + 0.5f) / (float)nZ;
          const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
          const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
          const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
          if (colRGBA[3] < 0.1f) {
            Solid[x][y][z]= true;
          }
          else {
            if (std::abs(colRGBA[0] - 0.5f) > 0.1f) PreBC[x][y][z]= true;
            if (std::abs(colRGBA[1] - 0.5f) > 0.1f) VelBC[x][y][z]= true;
            if (std::abs(colRGBA[2] - 0.5f) > 0.1f) SmoBC[x][y][z]= true;
          }
          if (PreBC[x][y][z]) {
            PresForced[x][y][z]= D.UI[BCPres__________].F() * ((colRGBA[0] > 0.5f) ? (1.0f) : (-1.0f));
          }
          if (VelBC[x][y][z]) {
            VelXForced[x][y][z]= D.UI[BCVelX__________].F() * ((colRGBA[1] > 0.5f) ? (1.0f) : (-1.0f));
            VelYForced[x][y][z]= D.UI[BCVelY__________].F() * ((colRGBA[1] > 0.5f) ? (1.0f) : (-1.0f));
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F() * ((colRGBA[1] > 0.5f) ? (1.0f) : (-1.0f));
          }
          if (SmoBC[x][y][z]) {
            SmokForced[x][y][z]= D.UI[BCSmok__________].F() * ((colRGBA[2] > 0.5f) ? (1.0f) : (-1.0f));
          }
        }
        // Double opposing inlets
        // |-----------|
        // |           |
        // | (>)   (<) |
        // |           |
        // |-----------|
        if (scenarioType == 1) {
          if (nY > 1 && (y == 0 || y == nY - 1)) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= 0.0f;
          }
          else if ((nX > 1 && (x == 0 || x == nX - 1)) ||
                   (nZ > 1 && (z == 0 || z == nZ - 1))) {
            Solid[x][y][z]= true;
          }
          else {
            const Vec::Vec3<float> posVox= Vec::Vec3<float>(D.boxMin[0], D.boxMin[1], D.boxMin[2]) + voxSize * Vec::Vec3<float>(x + 0.5f, y + 0.5f, z + 0.5f);
            const Vec::Vec3<float> posInlet(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
            for (int k= 0; k < 2; k++) {
              if (k == 0 && (posVox - posInlet).norm() > D.UI[ObjectSize0_____].F()) continue;
              if (k == 1 && (posVox - Vec::Vec3<float>(1.0f, 1.0f, 1.0f) + posInlet).norm() > D.UI[ObjectSize1_____].F()) continue;
              VelBC[x][y][z]= true;
              SmoBC[x][y][z]= true;
              VelXForced[x][y][z]= D.UI[BCVelX__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
              VelYForced[x][y][z]= D.UI[BCVelY__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
              VelZForced[x][y][z]= D.UI[BCVelZ__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
              SmokForced[x][y][z]= D.UI[BCSmok__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
            }
          }
        }
        // Circular obstacle in corridor showing vortex shedding
        // Test calib flow separation past cylinder https://link.springer.com/article/10.1007/s00521-020-05079-z
        // ---------------
        // >   O         >
        // ---------------
        if (scenarioType == 2) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
          }
          else if (y == nY - 1) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= 0.0f;
          }
          else if (y == 0) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
          else {
            const Vec::Vec3<float> posCell(((float)x + 0.5f) / (float)nX, ((float)y + 0.5f) / (float)nY, ((float)z + 0.5f) / (float)nZ);
            const Vec::Vec3<float> posObstacle(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
            Vec::Vec3<float> dist= (posCell - posObstacle);
            dist[0]*= (float)(nX - 1) * voxSize;
            dist[1]*= (float)(nY - 1) * voxSize;
            dist[2]*= (float)(nZ - 1) * voxSize;
            if (dist.norm() <= std::max(D.UI[ObjectSize0_____].F(), 0.0f))
              Solid[x][y][z]= true;
          }
        }
        // Cavity lid shear benchmark
        //  >>>>>>>>>>>
        // |           |
        // |           |
        // |           |
        // |-----------|
        if (scenarioType == 3) {
          if (y == 0 || y == nY - 1 || z == 0) {
            Solid[x][y][z]= true;
          }
          else if (z == nZ - 1) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
          }
          else if (y == nY / 2 && z > nZ / 2) {
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
        }
        // Flow constriction test with circular hole in a wall
        // -----|||--------
        // >    |||       >
        // >              >
        // >    |||       >
        // -----|||--------
        if (scenarioType == 4) {
          const Vec::Vec3<float> posVox= Vec::Vec3<float>(D.boxMin[0], D.boxMin[1], D.boxMin[2]) + voxSize * Vec::Vec3<float>(x + 0.5f, y + 0.5f, z + 0.5f);
          const Vec::Vec3<float> posWall(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
          const float radHole= D.UI[ObjectSize0_____].F();
          const float radWall= D.UI[ObjectSize1_____].F();
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            Solid[x][y][z]= true;
          }
          else if ((posVox - posWall).cwiseMul(Vec::Vec3<float>(0.0f, 1.0f, 0.0f)).norm() <= radWall &&
                   (posVox - posWall).cwiseMul(Vec::Vec3<float>(1.0f, 0.0f, 1.0f)).norm() > radHole) {
            Solid[x][y][z]= true;
          }
          else if (y == 0) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
          else if (y == nY - 1) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= 0.0f;
          }
        }
        // Central bloc with initial velocity
        // --------------
        // |    >>>>    |
        // |    >>>>    |
        // --------------
        if (scenarioType == 5) {
          if (((nX == 1) != (std::min(x, nX - 1 - x) > nX / 3)) &&
              ((nY == 1) != (std::min(y, nY - 1 - y) > nY / 3)) &&
              ((nZ == 1) != (std::min(z, nZ - 1 - z) > nZ / 3))) {
            VelX[x][y][z]= D.UI[BCVelX__________].F();
            VelY[x][y][z]= D.UI[BCVelY__________].F();
            VelZ[x][y][z]= D.UI[BCVelZ__________].F();
            Smok[x][y][z]= D.UI[BCSmok__________].F() * ((std::min(z, nZ - 1 - z) < 4 * (nZ - 1) / 9) ? (1.0f) : (-1.0f));
          }
        }
        // Poiseuille/Couette flow in tube with pressure gradient
        // ---------------
        // high        low
        // ---------------
        if (scenarioType == 6) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F() * ((z < nZ / 2) ? (-1.0f) : (1.0f));
            VelYForced[x][y][z]= D.UI[BCVelY__________].F() * ((z < nZ / 2) ? (-1.0f) : (1.0f));
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F() * ((z < nZ / 2) ? (-1.0f) : (1.0f));
          }
          else if (nY > 1 && (y == 0 || y == nY - 1)) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= D.UI[BCPres__________].F() * ((y < nY / 2) ? (1.0f) : (-1.0f));
          }
          else if (std::max(y, nY - 1 - y) == nY / 2) {
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
        }
        // Thermal convection cell
        // |---cold---|
        // |          |
        // |---warm---|
        if (scenarioType == 7) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nY > 1 && (y == 0 || y == nY - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            Solid[x][y][z]= true;
          }
          else if (nZ > 1 && std::min(z, nZ - 1 - z) < nZ / 3) {
            Smok[x][y][z]= D.UI[BCSmok__________].F() * ((z > nZ / 2) ? (-1.0f) : (1.0f));
            Smok[x][y][z]+= Random::Val(-0.01f, 0.01f);
          }
        }
        // Pipe of constant diameter with right angle turn and UI parameters to tweak position, curvature, diameters
        // ----------_
        // >          -
        // -----_      |
        //       |     |
        //       |  v  |
        if (scenarioType == 8) {
          const Vec::Vec3<float> posVox= Vec::Vec3<float>(D.boxMin[0], D.boxMin[1], D.boxMin[2]) + voxSize * Vec::Vec3<float>(x + 0.5f, y + 0.5f, z + 0.5f);
          const Vec::Vec3<float> posBend(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
          const float radPipe= D.UI[ObjectSize0_____].F();
          const float radBend= D.UI[ObjectSize1_____].F();
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nY > 1 && y == nY - 1) ||
              (nZ > 1 && z == nZ - 1)) {
            Solid[x][y][z]= true;
          }
          else if (posVox[1] < posBend[1]) {
            if ((posVox - posBend - Vec::Vec3<float>(0.0f, 0.0f, radBend + radPipe)).cwiseMul(Vec::Vec3<float>(1.0f, 0.0f, 1.0f)).norm() > radPipe) {
              Solid[x][y][z]= true;
            }
            else if (y == 0) {
              VelBC[x][y][z]= true;
              VelXForced[x][y][z]= D.UI[BCVelX__________].F();
              VelYForced[x][y][z]= D.UI[BCVelY__________].F();
              VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
              // PreBC[x][y][z]= true;
              // PresForced[x][y][z]= D.UI[BCPres__________].F();
              SmoBC[x][y][z]= true;
              SmokForced[x][y][z]= D.UI[BCSmok__________].F();
            }
          }
          else if (posVox[2] < posBend[2]) {
            if ((posVox - posBend - Vec::Vec3<float>(0.0f, radBend + radPipe, 0.0f)).cwiseMul(Vec::Vec3<float>(1.0f, 1.0f, 0.0f)).norm() > radPipe) {
              Solid[x][y][z]= true;
            }
            else if (z == 0) {
              PreBC[x][y][z]= true;
              PresForced[x][y][z]= 0.0;
            }
          }
          else if ((posBend + ((posVox - posBend).cwiseMul(Vec::Vec3<float>(0.0f, 1.0f, 1.0f)).normalized() * (radBend + radPipe)) - posVox).norm() > radPipe) {
            Solid[x][y][z]= true;
          }
        }
      }
    }
  }
}
