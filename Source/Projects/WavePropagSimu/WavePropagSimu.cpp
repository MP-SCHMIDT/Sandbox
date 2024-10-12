#include "WavePropagSimu.hpp"


// Standard lib
#include <array>
#include <format>
#include <numbers>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Draw/DrawField.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/BoxGrid.hpp"
#include "Type/Field.hpp"
#include "Math/Functions.hpp"
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"


// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
WavePropagSimu::WavePropagSimu() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void WavePropagSimu::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ScenarioPreset__", 0));
    D.UI.push_back(ParamUI("ScenarioFileID__", 0));
    D.UI.push_back(ParamUI("ResolutionX_____", 1));
    D.UI.push_back(ParamUI("ResolutionY_____", 200));
    D.UI.push_back(ParamUI("ResolutionZ_____", 200));
    D.UI.push_back(ParamUI("VoxelSize_______", 0.005));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("TimeStep________", 0.02));
    D.UI.push_back(ParamUI("Parallelize_____", 0));
    D.UI.push_back(ParamUI("NbSubsteps______", 1));
    D.UI.push_back(ParamUI("MaxWaveSpeed____", 0.05));
    D.UI.push_back(ParamUI("MaxAmplitude____", 1.0));
    D.UI.push_back(ParamUI("SourceRadius____", 0.05));
    D.UI.push_back(ParamUI("SourceFrequency_", 1.0));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("SliceDim________", 0));
    D.UI.push_back(ParamUI("SliceRelPosX____", 0.5));
    D.UI.push_back(ParamUI("SliceRelPosY____", 0.5));
    D.UI.push_back(ParamUI("SliceRelPosZ____", 0.5));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("ColorFactor_____", 8.0));
    D.UI.push_back(ParamUI("ScaleFactor_____", 0.1));
    D.UI.push_back(ParamUI("AlphaFactor_____", 0.0));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("TestParamWAV_0__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_1__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_2__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_3__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_4__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_5__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_6__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_7__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_8__", 0.0));
    D.UI.push_back(ParamUI("TestParamWAV_9__", 0.0));
    D.UI.push_back(ParamUI("______________03", NAN));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));

    D.displayModeLabel[1]= "Val";
    D.displayModeLabel[2]= "Densi";
    D.displayModeLabel[3]= "Elev";
    D.displayMode[1]= false;
    D.displayMode[2]= false;
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool WavePropagSimu::CheckAlloc() {
  if (D.UI[ScenarioPreset__].hasChanged()) isAllocated= false;
  if (D.UI[ScenarioFileID__].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionX_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionY_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionZ_____].hasChanged()) isAllocated= false;
  if (D.UI[VoxelSize_______].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool WavePropagSimu::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void WavePropagSimu::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  nX= std::max(D.UI[ResolutionX_____].I(), 1);
  nY= std::max(D.UI[ResolutionY_____].I(), 1);
  nZ= std::max(D.UI[ResolutionZ_____].I(), 1);
  simTime= 0;
  sourcePos= {NAN, NAN, NAN};

  UNew= Field::Field3(nX, nY, nZ, 0.0);
  UCur= Field::Field3(nX, nY, nZ, 0.0);
  UOld= Field::Field3(nX, nY, nZ, 0.0);
  Speed= Field::Field3(nX, nY, nZ, 1.0);
  BorderNeighbor= Field::Field3<std::array<char, 3>>(nX, nY, nZ, std::array<char, 3>{0, 0, 0});
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (nX >= 3 && x == 0) BorderNeighbor.at(x, y, z)[0]= +1;
        if (nY >= 3 && y == 0) BorderNeighbor.at(x, y, z)[1]= +1;
        if (nZ >= 3 && z == 0) BorderNeighbor.at(x, y, z)[2]= +1;
        if (nX >= 3 && x == nX - 1) BorderNeighbor.at(x, y, z)[0]= -1;
        if (nY >= 3 && y == nY - 1) BorderNeighbor.at(x, y, z)[1]= -1;
        if (nZ >= 3 && z == nZ - 1) BorderNeighbor.at(x, y, z)[2]= -1;
      }
    }
  }


  D.boxMin= {0.5 - 0.5 * (double)nX * D.UI[VoxelSize_______].D(),
             0.5 - 0.5 * (double)nY * D.UI[VoxelSize_______].D(),
             0.5 - 0.5 * (double)nZ * D.UI[VoxelSize_______].D()};
  D.boxMax= {0.5 + 0.5 * (double)nX * D.UI[VoxelSize_______].D(),
             0.5 + 0.5 * (double)nY * D.UI[VoxelSize_______].D(),
             0.5 + 0.5 * (double)nZ * D.UI[VoxelSize_______].D()};
}


// Refresh the project
void WavePropagSimu::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;


  // Get scenario ID and optionally load bitmap file
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (D.UI[ScenarioPreset__].I() == 0) {
    if (D.UI[ScenarioFileID__].I() == 0) FileInput::LoadImageBMPFile("./FileInput/Images/DoubleSlit.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 1) FileInput::LoadImageBMPFile("./FileInput/Images/FiberOptic.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 2) FileInput::LoadImageBMPFile("./FileInput/Images/LensConvex.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 3) FileInput::LoadImageBMPFile("./FileInput/Images/LensConcave.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 4) FileInput::LoadImageBMPFile("./FileInput/Images/Bimaterial.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 5) FileInput::LoadImageBMPFile("./FileInput/Images/SimpleSmile.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 6) FileInput::LoadImageBMPFile("./FileInput/Images/Logo.bmp", imageRGBA, false);
  }
  // Initialize scenario values
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Reset values
        UNew.at(x, y, z)= UCur.at(x, y, z)= UOld.at(x, y, z)= 0.0;
        Speed.at(x, y, z)= 1.0;

        // Scenario from loaded BMP file
        if (D.UI[ScenarioPreset__].I() == 0 && !imageRGBA.empty()) {
          const double posW= (double)(imageRGBA.size() - 1) * ((double)y + 0.5) / (double)nY;
          const double posH= (double)(imageRGBA[0].size() - 1) * ((double)z + 0.5) / (double)nZ;
          const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
          const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
          const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
          Speed.at(x, y, z)= (colRGBA[0] + colRGBA[1] + colRGBA[2]) / 3.0;
        }

        // Empty box domain with Dirichlet boundary
        if (D.UI[ScenarioPreset__].I() == 1) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nY > 1 && (y == 0 || y == nY - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1)))
            Speed.at(x, y, z)= 0.0;
        }

        // Simple sphere domain
        if (D.UI[ScenarioPreset__].I() == 2) {
          const int radius= std::min((nX > 1) ? nX : INT_MAX, std::min((nY > 1) ? nY : INT_MAX, (nZ > 1) ? nZ : INT_MAX)) / 2;
          if (std::pow(x - nX / 2, 2) + std::pow(y - nY / 2, 2) + std::pow(z - nZ / 2, 2) >= std::pow(radius - 2, 2)) Speed.at(x, y, z)= 0.0;
        }

        // Simple ellipse domain
        if (D.UI[ScenarioPreset__].I() == 3) {
          const int radius= std::min((nX > 1) ? nX : INT_MAX, std::min((nY > 1) ? nY : INT_MAX, (nZ > 1) ? nZ : INT_MAX)) / 2;
          if (2.0 * std::pow(x - nX / 2, 2) + std::pow(y - nY / 2, 2) + 2.0 * std::pow(z - nZ / 2, 2) >= std::pow(radius - 2, 2)) Speed.at(x, y, z)= 0.0;
        }

        // Convex lens domain as the intersection of two spheres
        if (D.UI[ScenarioPreset__].I() == 4) {
          const int radius= std::min((nX > 1) ? nX : INT_MAX, std::min((nY > 1) ? nY : INT_MAX, (nZ > 1) ? nZ : INT_MAX)) / 2;
          if (std::pow(x - nX / 2, 2) + std::pow(y - 1 * nY / 4, 2) + std::pow(z - nZ / 2, 2) < std::pow(radius - 2, 2) &&
              std::pow(x - nX / 2, 2) + std::pow(y - 3 * nY / 4, 2) + std::pow(z - nZ / 2, 2) < std::pow(radius - 2, 2)) Speed.at(x, y, z)= 0.5;
        }
      }
    }
  }
}


// Handle UI parameter change
void WavePropagSimu::ParamChange() {
}


// Handle keypress
void WavePropagSimu::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  if (D.keyLetterUpperCase == 'P') sourcePos= D.mouseProjX;
}


// Handle mouse action
void WavePropagSimu::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  if (D.mouseMiddleButtonState == 1 || D.mouseMiddleButtonState == 2) sourcePos= D.mouseProjX;
  if (D.mouseMiddleButtonState == 3) sourcePos= {NAN, NAN, NAN};
}


// Animate the project
void WavePropagSimu::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Simulate
  for (int k= 0; k < std::max(D.UI[NbSubsteps______].I(), 1); k++) {
    // Add source
    AddSource();

    // Update field values
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    StepSimulation();
    if (D.UI[VerboseLevel____].I() >= 1) printf("stepT %f\n", Timer::PopTimer());

    // Advance time
    simTime+= D.UI[TimeStep________].D();
  }

  // Draw the plot data
  D.Plot.resize(2);
  D.Plot[0].name= "USum";
  D.Plot[1].name= "USumAbs";
  for (int k= 0; k < (int)D.Plot.size(); k++)
    D.Plot[k].val.push_back(0.0);
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        D.Plot[0].val[D.Plot[0].val.size() - 1]+= UNew.at(x, y, z);
        D.Plot[1].val[D.Plot[1].val.size() - 1]+= std::abs(UNew.at(x, y, z));
      }
    }
  }

  // Write the status
  D.Status.clear();
  D.Status.resize(2);
  D.Status[0]= std::format("SimTime:{:.6f}s", simTime);
  D.Status[1]= std::format("CFL:{:<.6f}", D.UI[MaxWaveSpeed____].D() * D.UI[TimeStep________].D() / D.UI[VoxelSize_______].D());
}


// Draw the project
void WavePropagSimu::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the wave field
  if (D.displayMode[1] || D.displayMode[2]) {
    Field::Field3<char> Show(nX, nY, nZ, 1);
    Field::Field3<std::array<float, 4>> Color(nX, nY, nZ, {0.0, 0.0, 0.0, 1.0});
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (D.UI[SliceDim________].I() == 1 && x != (int)std::round(D.UI[SliceRelPosX____].D() * nX)) Show.at(x, y, z)= false;
          if (D.UI[SliceDim________].I() == 2 && y != (int)std::round(D.UI[SliceRelPosY____].D() * nY)) Show.at(x, y, z)= false;
          if (D.UI[SliceDim________].I() == 3 && z != (int)std::round(D.UI[SliceRelPosZ____].D() * nZ)) Show.at(x, y, z)= false;
          if (!Show.at(x, y, z)) continue;

          // Color by wave value
          if (D.displayMode[1]) {
            float r= 0.0, g= 0.0, b= 0.0;
            if (D.UI[ColorMode_______].I() == 0) Colormap::RatioToJetBrightSmooth(0.5 + UNew.at(x, y, z) * D.UI[ColorFactor_____].D(), r, g, b);
            if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToPlasma(std::abs(UNew.at(x, y, z)) * D.UI[ColorFactor_____].D(), r, g, b);
            Color.at(x, y, z)= {r, g, b, (float)std::abs(UNew.at(x, y, z))};
          }

          // Add shading to visualize wave speed field
          if (D.displayMode[2]) {
            float r= 0.0, g= 0.0, b= 0.0;
            Colormap::RatioToGrayscale(Speed.at(x, y, z), r, g, b);
            if (D.displayMode[1]) Color.at(x, y, z)= {0.6f * Color.at(x, y, z)[0] + 0.4f * r,
                                                      0.6f * Color.at(x, y, z)[1] + 0.4f * g,
                                                      0.6f * Color.at(x, y, z)[2] + 0.4f * b, Color.at(x, y, z)[3] + 0.1f * (1.0f - (float)Speed.at(x, y, z))};
            else Color.at(x, y, z)= {r, g, b, 0.1f * (1.0f - (float)Speed.at(x, y, z))};
          }

          // Apply transparency settings
          if (D.UI[AlphaFactor_____].F() > 0.0f)
            Color.at(x, y, z)[3]*= D.UI[AlphaFactor_____].F();
          else
            Color.at(x, y, z)[3]= 1.0f;
        }
      }
    }
    DrawField::DrawColored3DField(Show, Color, D.boxMin, D.boxMax, D.camDir, false, false, true);
  }

  // Draw the wave field as a 2.5D elevation map
  if (D.displayMode[3]) {
    if (nX == 1 && nY > 1 && nZ > 1) {
      const int x= 0;
      Field::Field2<Vec::Vec3<float>> terrainPos(nY, nZ, Vec::Vec3<float>(0.0, 0.0, 0.0));
      Field::Field2<Vec::Vec3<float>> terrainCol(nY, nZ, Vec::Vec3<float>(0.0, 0.0, 0.0));
      Field::Field2<Vec::Vec3<float>> terrainNor(nY, nZ, Vec::Vec3<float>(0.0, 0.0, 0.0));
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          terrainPos.at(y, z).set(0.5f * (D.boxMin[0] + D.boxMax[0]) + UNew.at(x, y, z) * D.UI[ScaleFactor_____].F(),
                                  D.boxMin[1] + ((float)y + 0.5f) * D.UI[VoxelSize_______].F(),
                                  D.boxMin[2] + ((float)z + 0.5f) * D.UI[VoxelSize_______].F());
          float r, g, b;
          if (D.UI[ColorMode_______].I() == 0) Colormap::RatioToViridis(0.5f + 0.5f * UNew.at(x, y, z) * D.UI[ColorFactor_____].F(), r, g, b);
          if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToPlasma(std::abs(UNew.at(x, y, z)) * D.UI[ColorFactor_____].F(), r, g, b);
          terrainCol.at(y, z).set(r, g, b);
        }
      }
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (y < nY - 1) {
            if (z < nZ - 1) terrainNor.at(y, z)+= ((terrainPos.at(y + 1, z) - terrainPos.at(y, z)).cross(terrainPos.at(y, z + 1) - terrainPos.at(y, z))).normalized();
            else if (z > 0) terrainNor.at(y, z)+= ((terrainPos.at(y, z - 1) - terrainPos.at(y, z)).cross(terrainPos.at(y + 1, z) - terrainPos.at(y, z))).normalized();
          }
          else if (y > 0) {
            if (z < nZ - 1) terrainNor.at(y, z)+= ((terrainPos.at(y, z + 1) - terrainPos.at(y, z)).cross(terrainPos.at(y - 1, z) - terrainPos.at(y, z))).normalized();
            else if (z > 0) terrainNor.at(y, z)+= ((terrainPos.at(y - 1, z) - terrainPos.at(y, z)).cross(terrainPos.at(y, z - 1) - terrainPos.at(y, z))).normalized();
          }
          terrainNor.at(y, z).normalize();
        }
      }
      glEnable(GL_LIGHTING);
      glBegin(GL_QUADS);
      for (int y= 0; y < nY - 1; y++) {
        for (int z= 0; z < nZ - 1; z++) {
          if (Speed.at(x, y, z) == 0.0 && Speed.at(x, y, z + 1) == 0.0 && Speed.at(x, y + 1, z) == 0.0 && Speed.at(x, y + 1, z + 1) == 0.0) continue;
          Vec::Vec3<float> flatNormal= (terrainNor.at(y, z) + terrainNor.at(y + 1, z) + terrainNor.at(y + 1, z + 1) + terrainNor.at(y, z + 1)).normalized();
          Vec::Vec3<float> flatColor= (terrainCol.at(y, z) + terrainCol.at(y + 1, z) + terrainCol.at(y + 1, z + 1) + terrainCol.at(y, z + 1)) / 4.0f;
          glColor3fv(flatColor.array());
          glNormal3fv(flatNormal.array());
          glVertex3fv(terrainPos.at(y, z).array());
          glVertex3fv(terrainPos.at(y + 1, z).array());
          glVertex3fv(terrainPos.at(y + 1, z + 1).array());
          glVertex3fv(terrainPos.at(y, z + 1).array());
        }
      }
      glEnd();
      glDisable(GL_LIGHTING);
    }
  }
}


void WavePropagSimu::StepSimulation() {
  UOld.swap(UCur);
  UCur.swap(UNew);
  // Explicit numerical integration of the wave equation
  // ∂²u/∂t² = c² ∂²u/∂x²
#pragma omp parallel for collapse(3) if (D.UI[Parallelize_____].I() > 0)
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Speed.at(x, y, z) == 0.0) continue;
        // Precompute data
        const double gamma= Speed.at(x, y, z) * D.UI[MaxWaveSpeed____].D() * D.UI[TimeStep________].D() / D.UI[VoxelSize_______].D();

        // Apply the absorbing boundary conditions with Perfectly Matched Layer on domain faces
        // https://en.wikipedia.org/wiki/Perfectly_matched_layer
        // https://hal.science/hal-01374183
        if (BorderNeighbor.at(x, y, z) != std::array<char, 3>{0, 0, 0}) {
          const std::array<char, 3> offset= BorderNeighbor.at(x, y, z);
          UNew.at(x, y, z)= UCur.at(x, y, z) - gamma * (UCur.at(x, y, z) - UCur.at(x + offset[0], y + offset[1], z + offset[2]));
        }
        else {
          // Compute the 3D field Laplacian with spatial discretization
          // https://en.wikipedia.org/wiki/Discrete_Laplace_operator
          double laplacian= 0.0;
          laplacian+= (x - 1 >= 0) ? UCur.at(x - 1, y, z) : UCur.at(x, y, z);
          laplacian+= (y - 1 >= 0) ? UCur.at(x, y - 1, z) : UCur.at(x, y, z);
          laplacian+= (z - 1 >= 0) ? UCur.at(x, y, z - 1) : UCur.at(x, y, z);
          laplacian+= (x + 1 < nX) ? UCur.at(x + 1, y, z) : UCur.at(x, y, z);
          laplacian+= (y + 1 < nY) ? UCur.at(x, y + 1, z) : UCur.at(x, y, z);
          laplacian+= (z + 1 < nZ) ? UCur.at(x, y, z + 1) : UCur.at(x, y, z);
          laplacian-= 6.0 * UCur.at(x, y, z);

          // Update the field values with time discretized wave equation
          // https://en.wikipedia.org/wiki/Wave_equation
          // https://vitalitylearning.medium.com/solving-the-1d-wave-equation-numerical-discretization-190a92c917bc
          // https://www.youtube.com/watch?v=pN-gi_omIVE
          // https://www.idpoisson.fr/berglund/wave_billiard.c
          UNew.at(x, y, z)= 2.0 * UCur.at(x, y, z) - UOld.at(x, y, z) + gamma * gamma * laplacian;
        }
      }
    }
  }
}


void WavePropagSimu::AddSource() {
  if (!std::isnan(sourcePos[0])) {
    const Vec::Vec3 source(sourcePos[0], sourcePos[1], sourcePos[2]);
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (Speed.at(x, y, z) == 0.0) continue;
          const Vec::Vec3 pos(D.boxMin[0] + ((double)x + 0.5) * D.UI[VoxelSize_______].D(),
                              D.boxMin[1] + ((double)y + 0.5) * D.UI[VoxelSize_______].D(),
                              D.boxMin[2] + ((double)z + 0.5) * D.UI[VoxelSize_______].D());
          const double blend= Functions::SmoothStep((D.UI[SourceRadius____].D() - (pos - source).norm()) / (2.0 * D.UI[SourceRadius____].D()));
          const double omega= D.UI[SourceFrequency_].D() * simTime * std::numbers::pi;
          UNew.at(x, y, z)= blend * std::sin(omega) * D.UI[MaxAmplitude____].D() + (1.0 - blend) * UNew.at(x, y, z);
          UCur.at(x, y, z)= blend * std::sin(omega) * D.UI[MaxAmplitude____].D() + (1.0 - blend) * UCur.at(x, y, z);
          UOld.at(x, y, z)= blend * std::sin(omega) * D.UI[MaxAmplitude____].D() + (1.0 - blend) * UOld.at(x, y, z);
        }
      }
    }
  }
}
