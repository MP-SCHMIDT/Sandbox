#include "FractalElevMap.hpp"


// Standard lib
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
FractalElevMap::FractalElevMap() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void FractalElevMap::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ResolutionX_____", 1000.0));       // Resolution of 2.5D height map
    D.UI.push_back(ParamUI("ResolutionY_____", 1000.0));       // Resolution of 2.5D height map
    D.UI.push_back(ParamUI("ZoomLevel_______", 0.5));          // Zoom factor
    D.UI.push_back(ParamUI("NumItera________", 50.0));         // Number of iterations of complex square + complex add
    D.UI.push_back(ParamUI("ShiftX__________", 0.365242330));  // Hand picked param for cool fractal zoom
    D.UI.push_back(ParamUI("ShiftY__________", 0.534752180));  // Hand picked param for cool fractal zoom
    D.UI.push_back(ParamUI("CoeffA__________", -0.8350));      // Hand picked param for cool fractal zoom
    D.UI.push_back(ParamUI("CoeffB__________", -0.2241));      // Hand picked param for cool fractal zoom
    D.UI.push_back(ParamUI("HeightLvls______", 32.0));         // Divergence threshold of complex iterations
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
bool FractalElevMap::CheckAlloc() {
  if (D.UI[ResolutionX_____].hasChanged()) isAllocated= false;
  if (D.UI[ResolutionY_____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool FractalElevMap::CheckRefresh() {
  if (D.UI[ZoomLevel_______].hasChanged()) isRefreshed= false;
  if (D.UI[NumItera________].hasChanged()) isRefreshed= false;
  if (D.UI[ShiftX__________].hasChanged()) isRefreshed= false;
  if (D.UI[ShiftY__________].hasChanged()) isRefreshed= false;
  if (D.UI[CoeffA__________].hasChanged()) isRefreshed= false;
  if (D.UI[CoeffB__________].hasChanged()) isRefreshed= false;
  if (D.UI[HeightLvls______].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void FractalElevMap::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  mapNbX= std::max(D.UI[ResolutionX_____].I(), 2);
  mapNbY= std::max(D.UI[ResolutionY_____].I(), 2);

  // Allocate data
  mapPos= Field::AllocField2D(mapNbX, mapNbY, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  mapNor= Field::AllocField2D(mapNbX, mapNbY, Vec::Vec3<float>(0.0f, 0.0f, 1.0f));
  mapCol= Field::AllocField2D(mapNbX, mapNbY, Vec::Vec3<float>(0.5f, 0.5f, 0.5f));
}


// Refresh the project
void FractalElevMap::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  mapZoom= std::max(D.UI[ZoomLevel_______].D(), 1.e-6);

  mapNbIter= std::max(D.UI[NumItera________].I(), 1);

  mapFocus= Vec::Vec2<double>(D.UI[ShiftX__________].D(), D.UI[ShiftY__________].D());
  mapConst= Vec::Vec2<double>(D.UI[CoeffA__________].D(), D.UI[CoeffB__________].D());

  mapDivThresh= std::max(D.UI[HeightLvls______].D(), 0.0);


  // Compute positions
#pragma omp parallel for
  for (int x= 0; x < mapNbX; x++) {
    for (int y= 0; y < mapNbY; y++) {
      mapPos[x][y][0]= float(x) / float(mapNbX - 1);
      mapPos[x][y][1]= float(y) / float(mapNbY - 1);

      Vec::Vec2<double> z= mapFocus + Vec::Vec2<double>(2.0 * double(x) / double(mapNbX - 1) - 1.0, 2.0 * double(y) / double(mapNbY - 1) - 1.0) / mapZoom;
      int idxIter= 0;
      while (idxIter < mapNbIter && z.normSquared() < mapDivThresh) {
        const double zr= z[0] * z[0] - z[1] * z[1];
        const double zi= 2.0 * z[0] * z[1];
        z= Vec::Vec2<double>(zr, zi) + mapConst;
        idxIter++;
      }

      // double val= double(idxIter) / double(mapNbIter-1);
      double val= (double(idxIter) - std::log2(std::max(std::log2(z.normSquared()), 1.0))) / double(mapNbIter - 1);
      // double val= - std::log2(std::max(std::log2(z.normSquared()), 1.0));
      // if (val != 0.0) val= std::log2(std::max(std::log2(val), 1.0));

      // mapPos[x][y][2]= float(z.norm());
      // if (mapPos[x][y][2] != mapPos[x][y][2]) mapPos[x][y][2]= D.UI[HeightLvls__].Get();
      Colormap::RatioToJetSmooth(float(val), mapCol[x][y][0], mapCol[x][y][1], mapCol[x][y][2]);

      mapPos[x][y][2]= 0.5f + 0.04f * std::min(std::max(float(val), 0.0f), 1.0f);
    }
  }

  // Smooth the positions
  for (int iter= 0; iter < std::max(mapNbX, mapNbY) / 128; iter++) {
    std::vector<std::vector<Vec::Vec3<float>>> mapPosOld= mapPos;
#pragma omp parallel for
    for (int x= 0; x < mapNbX; x++) {
      for (int y= 0; y < mapNbY; y++) {
        int count= 0;
        mapPos[x][y][2]= 0.0;
        for (int xOff= std::max(x - 1, 0); xOff <= std::min(x + 1, mapNbX - 1); xOff++) {
          for (int yOff= std::max(y - 1, 0); yOff <= std::min(y + 1, mapNbY - 1); yOff++) {
            mapPos[x][y][2]+= mapPosOld[xOff][yOff][2];
            count++;
          }
        }
        mapPos[x][y][2]/= float(count);
      }
    }
  }

  // Compute normals
#pragma omp parallel for
  for (int x= 0; x < mapNbX; x++) {
    for (int y= 0; y < mapNbY; y++) {
      mapNor[x][y].set(0.0f, 0.0f, 0.0f);
      if (x > 0 && y > 0)
        mapNor[x][y]+= ((mapPos[x - 1][y] - mapPos[x][y]).cross(mapPos[x][y - 1] - mapPos[x][y])).normalized();
      if (x < mapNbX - 1 && y > 0)
        mapNor[x][y]+= ((mapPos[x][y - 1] - mapPos[x][y]).cross(mapPos[x + 1][y] - mapPos[x][y])).normalized();
      if (x < mapNbX - 1 && y < mapNbY - 1)
        mapNor[x][y]+= ((mapPos[x + 1][y] - mapPos[x][y]).cross(mapPos[x][y + 1] - mapPos[x][y])).normalized();
      if (x > 0 && y < mapNbY - 1)
        mapNor[x][y]+= ((mapPos[x][y + 1] - mapPos[x][y]).cross(mapPos[x - 1][y] - mapPos[x][y])).normalized();
      mapNor[x][y].normalize();
    }
  }
}


// Handle keypress
void FractalElevMap::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Animate the project
void FractalElevMap::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
}


// Draw the project
void FractalElevMap::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the map
  glEnable(GL_LIGHTING);
  glBegin(GL_QUADS);
  for (int x= 0; x < mapNbX - 1; x++) {
    for (int y= 0; y < mapNbY - 1; y++) {
      Vec::Vec3<float> flatNormal= (mapNor[x][y] + mapNor[x + 1][y] + mapNor[x + 1][y + 1] + mapNor[x][y + 1]).normalized();
      Vec::Vec3<float> flatColor= (mapCol[x][y] + mapCol[x + 1][y] + mapCol[x + 1][y + 1] + mapCol[x][y + 1]) / 4.0f;
      glColor3fv(flatColor.array());
      glNormal3fv(flatNormal.array());
      glVertex3fv(mapPos[x][y].array());
      glVertex3fv(mapPos[x + 1][y].array());
      glVertex3fv(mapPos[x + 1][y + 1].array());
      glVertex3fv(mapPos[x][y + 1].array());
    }
  }
  glEnd();
  glDisable(GL_LIGHTING);
}
