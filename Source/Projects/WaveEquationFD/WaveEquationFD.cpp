#include "WaveEquationFD.hpp"


// Standard lib
#include <array>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Geom/BoxGrid.hpp"
#include "Math/Field.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
WaveEquationFD::WaveEquationFD() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void WaveEquationFD::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ResolutionX_____", 1));
    D.UI.push_back(ParamUI("ResolutionY_____", 100));
    D.UI.push_back(ParamUI("ResolutionZ_____", 100));
    D.UI.push_back(ParamUI("VoxelSize_______", 0.01));
    D.UI.push_back(ParamUI("TimeStep________", 0.01));
    D.UI.push_back(ParamUI("MaxWaveSpeed____", 0.01));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("Multithread_____", 0));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("______________01", NAN));
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
    D.UI.push_back(ParamUI("______________02", NAN));
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
bool WaveEquationFD::CheckAlloc() {
  if (D.UI[ResolutionX_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionY_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionZ_____].hasChanged()) isAllocated= false;
  if (D.UI[VoxelSize_______].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool WaveEquationFD::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void WaveEquationFD::Allocate() {
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
void WaveEquationFD::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Initialize scenario values
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        UNew[x][y][z]= 0.0;
        UCur[x][y][z]= 0.0;
        UOld[x][y][z]= 0.0;
        Speed[x][y][z]= 1.0;
        if (nX > 1 && (x == 0 || x == nX - 1)) Speed[x][y][z]= 0.0;
        if (nY > 1 && (y == 0 || y == nY - 1)) Speed[x][y][z]= 0.0;
        if (nZ > 1 && (z == 0 || z == nZ - 1)) Speed[x][y][z]= 0.0;
      }
    }
  }
}


// Handle keypress
void WaveEquationFD::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Handle mouse action
void WaveEquationFD::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable

  // Get cursor coordinates in the field
  const int xCursor= (int)std::floor((double)nX * (D.mouseProjX[0] - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
  const int yCursor= (int)std::floor((double)nY * (D.mouseProjX[1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
  const int zCursor= (int)std::floor((double)nZ * (D.mouseProjX[2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
  if (xCursor < 0 || xCursor >= nX) return;
  if (yCursor < 0 || yCursor >= nY) return;
  if (zCursor < 0 || zCursor >= nZ) return;

  // Set the target voxel
  if (D.UI[TestParamWAV_0__].B()) UNew[xCursor][yCursor][zCursor]= 1.0;
  if (D.UI[TestParamWAV_1__].B()) UCur[xCursor][yCursor][zCursor]= 1.0;
  if (D.UI[TestParamWAV_2__].B()) UOld[xCursor][yCursor][zCursor]= 1.0;
}


// Animate the project
void WaveEquationFD::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Update field values
  UOld= UCur;
  UCur= UNew;
#pragma omp parallel for collapse(3) if (D.UI[Multithread_____].B())
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        UNew[x][y][z]= 0.0;
        if (Speed[x][y][z] == 0.0) continue;
        int count= 0;
        double sum= 0.0;
        if (x - 1 >= 0 && ++count) sum+= UCur[x - 1][y][z];
        if (y - 1 >= 0 && ++count) sum+= UCur[x][y - 1][z];
        if (z - 1 >= 0 && ++count) sum+= UCur[x][y][z - 1];
        if (x + 1 < nX && ++count) sum+= UCur[x + 1][y][z];
        if (y + 1 < nY && ++count) sum+= UCur[x][y + 1][z];
        if (z + 1 < nZ && ++count) sum+= UCur[x][y][z + 1];
        const double gamma= D.UI[MaxWaveSpeed____].D() * Speed[x][y][z] * D.UI[TimeStep________].D() / D.UI[VoxelSize_______].D();
        UNew[x][y][z]= 2.0 * UCur[x][y][z] - UOld[x][y][z] + gamma * gamma * (sum - (double)count * UCur[x][y][z]);
      }
    }
  }

  // Advance time
  simTime+= D.UI[TimeStep________].D();
}


// Draw the project
void WaveEquationFD::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;


  // Draw the scalar fields
  if (D.displayMode2) {
    // Set the scene transformation
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5 * D.UI[VoxelSize_______].D(),
                 D.boxMin[1] + 0.5 * D.UI[VoxelSize_______].D(),
                 D.boxMin[2] + 0.5 * D.UI[VoxelSize_______].D());
    glScalef(D.UI[VoxelSize_______].D(), D.UI[VoxelSize_______].D(), D.UI[VoxelSize_______].D());
    if (nX == 1) glScalef(0.1, 1.0, 1.0);
    if (nY == 1) glScalef(1.0, 0.1, 1.0);
    if (nZ == 1) glScalef(1.0, 1.0, 0.1);
    // Sweep the field
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          float r= 0.0, g= 0.0, b= 0.0;
          Colormap::RatioToViridis(0.5 + UCur[x][y][z] * D.UI[ColorFactor_____].D(), r, g, b);
          glColor3f(r, g, b);
          glPushMatrix();
          glTranslatef((double)x, (double)y, (double)z);
          glutSolidCube(1.0);
          glPopMatrix();
        }
      }
    }
    glPopMatrix();
  }
}
