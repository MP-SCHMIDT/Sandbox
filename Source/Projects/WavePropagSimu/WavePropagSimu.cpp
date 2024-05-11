#include "WavePropagSimu.hpp"


// Standard lib
#include <array>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Draw/DrawField.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/BoxGrid.hpp"
#include "Math/Field.hpp"
#include "Math/Functions.hpp"
#include "Math/Vec.hpp"


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
    D.UI.push_back(ParamUI("ResolutionY_____", 100));
    D.UI.push_back(ParamUI("ResolutionZ_____", 100));
    D.UI.push_back(ParamUI("VoxelSize_______", 0.01));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("TimeStep________", 0.05));
    D.UI.push_back(ParamUI("Parallelize_____", 0));
    D.UI.push_back(ParamUI("NbSubsteps______", 1));
    D.UI.push_back(ParamUI("MaxWaveSpeed____", 0.05));
    D.UI.push_back(ParamUI("MaxAmplitude____", 1.0));
    D.UI.push_back(ParamUI("BrushRadius_____", 0.05));
    D.UI.push_back(ParamUI("BrushBorder_____", 0.04));
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

  UNew= Field::AllocField3D(nX, nY, nZ, 0.0);
  UCur= Field::AllocField3D(nX, nY, nZ, 0.0);
  UOld= Field::AllocField3D(nX, nY, nZ, 0.0);
  Speed= Field::AllocField3D(nX, nY, nZ, 1.0);

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
    if (D.UI[ScenarioFileID__].I() == 1) FileInput::LoadImageBMPFile("./FileInput/Images/LensConvex.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 2) FileInput::LoadImageBMPFile("./FileInput/Images/LensConcave.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 3) FileInput::LoadImageBMPFile("./FileInput/Images/Bimaterial.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 4) FileInput::LoadImageBMPFile("./FileInput/Images/SimpleSmile.bmp", imageRGBA, false);
    if (D.UI[ScenarioFileID__].I() == 5) FileInput::LoadImageBMPFile("./FileInput/Images/Logo.bmp", imageRGBA, false);
  }
  // Initialize scenario values
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Reset values
        UNew[x][y][z]= 0.0;
        UCur[x][y][z]= 0.0;
        UOld[x][y][z]= 0.0;
        Speed[x][y][z]= 1.0;

        // Scenario from loaded BMP file
        if (D.UI[ScenarioPreset__].I() == 0 && !imageRGBA.empty()) {
          const double posW= (double)(imageRGBA.size() - 1) * ((double)y + 0.5) / (double)nY;
          const double posH= (double)(imageRGBA[0].size() - 1) * ((double)z + 0.5) / (double)nZ;
          const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
          const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
          const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
          Speed[x][y][z]= (colRGBA[0] + colRGBA[1] + colRGBA[2]) / 3.0;
        }

        // Empty box domain with Dirichlet boundary
        if (D.UI[ScenarioPreset__].I() == 1) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nY > 1 && (y == 0 || y == nY - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1)))
            Speed[x][y][z]= 0.0;
        }

        // Simple sphere domain
        if (D.UI[ScenarioPreset__].I() == 2) {
          const int radius= std::max(nX, std::max(nY, nZ)) / 2;
          if (std::pow(x - nX / 2, 2) + std::pow(y - nY / 2, 2) + std::pow(z - nZ / 2, 2) >= std::pow(radius - 2, 2)) Speed[x][y][z]= 0.0;
        }

        // Simple ellipse domain
        if (D.UI[ScenarioPreset__].I() == 3) {
          const int radius= std::max(nX, std::max(nY, nZ)) / 2;
          if (2.0 * std::pow(x - nX / 2, 2) + std::pow(y - nY / 2, 2) + 2.0 * std::pow(z - nZ / 2, 2) >= std::pow(radius - 2, 2)) Speed[x][y][z]= 0.0;
        }
      }
    }
  }
}


// Handle keypress
void WavePropagSimu::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Handle mouse action
void WavePropagSimu::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable

  // Set the target voxels
  const Vec::Vec3 cursor(D.mouseProjX[0], D.mouseProjX[1], D.mouseProjX[2]);
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        const Vec::Vec3 pos(D.boxMin[0] + (double(x) + 0.5) * D.UI[VoxelSize_______].D(),
                            D.boxMin[1] + (double(y) + 0.5) * D.UI[VoxelSize_______].D(),
                            D.boxMin[2] + (double(z) + 0.5) * D.UI[VoxelSize_______].D());
        const double dist= D.UI[BrushRadius_____].D() - (pos - cursor).norm();
        const double val= 2.0 * Functions::SmoothHeaviside(dist, D.UI[BrushBorder_____].D()) - 1.0;
        UNew[x][y][z]= std::max(D.UI[MaxAmplitude____].D() * val, UNew[x][y][z]);
        UCur[x][y][z]= std::max(D.UI[MaxAmplitude____].D() * val, UCur[x][y][z]);
        UOld[x][y][z]= std::max(D.UI[MaxAmplitude____].D() * val, UOld[x][y][z]);
      }
    }
  }
}


// Animate the project
void WavePropagSimu::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  for (int k= 0; k < std::max(D.UI[NbSubsteps______].I(), 1); k++) {
    // Update field values
    StepSimulation();

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
        D.Plot[0].val[D.Plot[0].val.size() - 1]+= UNew[x][y][z];
        D.Plot[1].val[D.Plot[1].val.size() - 1]+= std::abs(UNew[x][y][z]);
      }
    }
  }

  // Write the status
  if ((int)D.Status.size() != 2) D.Status.resize(2);
  D.Status[0]= std::string{"SimTime: "} + std::to_string(simTime);
  D.Status[1]= std::string{"CFL: "} + std::to_string(D.UI[MaxWaveSpeed____].D() * D.UI[TimeStep________].D() / D.UI[VoxelSize_______].D());
}


// Draw the project
void WavePropagSimu::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the wave field
  if (D.displayMode1 || D.displayMode2) {
    std::vector<std::vector<std::vector<bool>>> Show= Field::AllocField3D(nX, nY, nZ, true);
    std::vector<std::vector<std::vector<std::array<float, 4>>>> Color= Field::AllocField3D(nX, nY, nZ, std::array<float, 4>{0.0, 0.0, 0.0, 1.0});
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (D.UI[SliceDim________].I() == 1 && x != (int)std::round(D.UI[SliceRelPosX____].D() * nX)) Show[x][y][z]= false;
          if (D.UI[SliceDim________].I() == 2 && y != (int)std::round(D.UI[SliceRelPosY____].D() * nY)) Show[x][y][z]= false;
          if (D.UI[SliceDim________].I() == 3 && z != (int)std::round(D.UI[SliceRelPosZ____].D() * nZ)) Show[x][y][z]= false;
          if (!Show[x][y][z]) continue;

          // Color by wave value
          if (D.displayMode1) {
            float r= 0.0, g= 0.0, b= 0.0;
            if (D.UI[ColorMode_______].I() == 0) Colormap::RatioToViridis(0.5 + UNew[x][y][z] * D.UI[ColorFactor_____].D(), r, g, b);
            if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToPlasma(std::abs(UNew[x][y][z]) * D.UI[ColorFactor_____].D(), r, g, b);
            Color[x][y][z]= {r, g, b, (float)std::abs(UNew[x][y][z])};
          }

          // Add shading to visualize wave speed field
          if (D.displayMode2) {
            float r= 0.0, g= 0.0, b= 0.0;
            Colormap::RatioToGrayscale(Speed[x][y][z], r, g, b);
            if (D.displayMode1) Color[x][y][z]= {0.6f * Color[x][y][z][0] + 0.4f * r,
                                                 0.6f * Color[x][y][z][1] + 0.4f * g,
                                                 0.6f * Color[x][y][z][2] + 0.4f * b, Color[x][y][z][3]};
            else Color[x][y][z]= {r, g, b, Color[x][y][z][2]};
          }

          // Apply transparency settings
          if (D.UI[AlphaFactor_____].F() > 0.0f)
            Color[x][y][z][3]*= D.UI[AlphaFactor_____].F();
          else
            Color[x][y][z][3]= 1.0f;
        }
      }
    }
    DrawField::DrawColored3DField(Show, Color, D.boxMin, D.boxMax, D.camDir, false, false, true);
  }

  // Draw the wave field as a 2.5D elevation map
  if (D.displayMode3) {
    if (nX == 1 && nY > 1 && nZ > 1) {
      std::vector<std::vector<Vec::Vec3<float>>> terrainPos= Field::AllocField2D(nY, nZ, Vec::Vec3<float>(0.0, 0.0, 0.0));
      std::vector<std::vector<Vec::Vec3<float>>> terrainCol= Field::AllocField2D(nY, nZ, Vec::Vec3<float>(0.0, 0.0, 0.0));
      std::vector<std::vector<Vec::Vec3<float>>> terrainNor= Field::AllocField2D(nY, nZ, Vec::Vec3<float>(0.0, 0.0, 0.0));
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          terrainPos[y][z].set(0.5f * (D.boxMin[0] + D.boxMax[0]) + UNew[0][y][z] * D.UI[ScaleFactor_____].F(),
                               D.boxMin[1] + ((float)y + 0.5f) * D.UI[VoxelSize_______].F(),
                               D.boxMin[2] + ((float)z + 0.5f) * D.UI[VoxelSize_______].F());
          float r, g, b;
          if (D.UI[ColorMode_______].I() == 0) Colormap::RatioToViridis(0.5f + 0.5f * UNew[0][y][z] * D.UI[ColorFactor_____].F(), r, g, b);
          if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToPlasma(std::abs(UNew[0][y][z]) * D.UI[ColorFactor_____].F(), r, g, b);
          terrainCol[y][z].set(r, g, b);
        }
      }
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (y < nY - 1) {
            if (z < nZ - 1) terrainNor[y][z]+= ((terrainPos[y + 1][z] - terrainPos[y][z]).cross(terrainPos[y][z + 1] - terrainPos[y][z])).normalized();
            else if (z > 0) terrainNor[y][z]+= ((terrainPos[y][z - 1] - terrainPos[y][z]).cross(terrainPos[y + 1][z] - terrainPos[y][z])).normalized();
          }
          else if (y > 0) {
            if (z < nZ - 1) terrainNor[y][z]+= ((terrainPos[y][z + 1] - terrainPos[y][z]).cross(terrainPos[y - 1][z] - terrainPos[y][z])).normalized();
            else if (z > 0) terrainNor[y][z]+= ((terrainPos[y - 1][z] - terrainPos[y][z]).cross(terrainPos[y][z - 1] - terrainPos[y][z])).normalized();
          }
          terrainNor[y][z].normalize();
        }
      }
      glEnable(GL_LIGHTING);
      glBegin(GL_QUADS);
      for (int y= 0; y < nY - 1; y++) {
        for (int z= 0; z < nZ - 1; z++) {
          if (Speed[0][y][z] == 0.0 && Speed[0][y][z + 1] == 0.0 && Speed[0][y + 1][z] == 0.0 && Speed[0][y + 1][z + 1] == 0.0) continue;
          Vec::Vec3<float> flatNormal= (terrainNor[y][z] + terrainNor[y + 1][z] + terrainNor[y + 1][z + 1] + terrainNor[y][z + 1]).normalized();
          Vec::Vec3<float> flatColor= (terrainCol[y][z] + terrainCol[y + 1][z] + terrainCol[y + 1][z + 1] + terrainCol[y][z + 1]) / 4.0f;
          glColor3fv(flatColor.array());
          glNormal3fv(flatNormal.array());
          glVertex3fv(terrainPos[y][z].array());
          glVertex3fv(terrainPos[y + 1][z].array());
          glVertex3fv(terrainPos[y + 1][z + 1].array());
          glVertex3fv(terrainPos[y][z + 1].array());
        }
      }
      glEnd();
      glDisable(GL_LIGHTING);
    }
  }
}


void WavePropagSimu::StepSimulation() {
  UOld= UCur;
  UCur= UNew;
#pragma omp parallel for collapse(3) if (D.UI[Parallelize_____].B())
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Reset and precompute
        UNew[x][y][z]= 0.0;
        if (Speed[x][y][z] == 0.0) continue;
        const double gamma= D.UI[MaxWaveSpeed____].D() * Speed[x][y][z] * D.UI[TimeStep________].D() / D.UI[VoxelSize_______].D();

        // Update the field values with dicretized wave equation
        double sum= 0.0;
        sum+= (x - 1 >= 0) ? (UCur[x - 1][y][z]) : (UCur[x][y][z]);
        sum+= (y - 1 >= 0) ? (UCur[x][y - 1][z]) : (UCur[x][y][z]);
        sum+= (z - 1 >= 0) ? (UCur[x][y][z - 1]) : (UCur[x][y][z]);
        sum+= (x + 1 < nX) ? (UCur[x + 1][y][z]) : (UCur[x][y][z]);
        sum+= (y + 1 < nY) ? (UCur[x][y + 1][z]) : (UCur[x][y][z]);
        sum+= (z + 1 < nZ) ? (UCur[x][y][z + 1]) : (UCur[x][y][z]);
        UNew[x][y][z]= 2.0 * UCur[x][y][z] - UOld[x][y][z] + gamma * gamma * (sum - 6.0 * UCur[x][y][z]);

        // Apply the absorbing boundary conditions with Perfectly Matched Layer on domain faces
        // https://en.wikipedia.org/wiki/Perfectly_matched_layer
        // https://hal.science/hal-01374183
        // https://www.idpoisson.fr/berglund/wave_billiard.c
        if (nX >= 3 && x == nX - 1) UNew[x][y][z]= UCur[x][y][z] - gamma * (UCur[x][y][z] - UCur[x - 1][y][z]);
        else if (nX >= 3 && x == 0) UNew[x][y][z]= UCur[x][y][z] - gamma * (UCur[x][y][z] - UCur[x + 1][y][z]);
        if (nY >= 3 && y == nY - 1) UNew[x][y][z]= UCur[x][y][z] - gamma * (UCur[x][y][z] - UCur[x][y - 1][z]);
        else if (nY >= 3 && y == 0) UNew[x][y][z]= UCur[x][y][z] - gamma * (UCur[x][y][z] - UCur[x][y + 1][z]);
        if (nZ >= 3 && z == nZ - 1) UNew[x][y][z]= UCur[x][y][z] - gamma * (UCur[x][y][z] - UCur[x][y][z - 1]);
        else if (nZ >= 3 && z == 0) UNew[x][y][z]= UCur[x][y][z] - gamma * (UCur[x][y][z] - UCur[x][y][z + 1]);
      }
    }
  }
}
