#include "TerrainErosion.hpp"


// Standard lib
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Type/Field.hpp"
#include "Type/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
TerrainErosion::TerrainErosion() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void TerrainErosion::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("TerrainNbX______", 128));
    D.UI.push_back(ParamUI("TerrainNbY______", 128));
    D.UI.push_back(ParamUI("TerrainNbCut____", 256));
    D.UI.push_back(ParamUI("DropletNbK______", 1000));
    D.UI.push_back(ParamUI("DropletRad______", 0.01));
    D.UI.push_back(ParamUI("SimuTimestep____", 0.02));
    D.UI.push_back(ParamUI("VelDecay________", 0.5));
    D.UI.push_back(ParamUI("ErosionCoeff____", 0.05));
    D.UI.push_back(ParamUI("SmoothResist____", 0.99));
    D.UI.push_back(ParamUI("CliffThresh_____", 0.80));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));

    D.displayModeLabel[1]= "Terrain";
    D.displayModeLabel[4]= "Normal";
    D.displayModeLabel[5]= "Drop";
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool TerrainErosion::CheckAlloc() {
  if (D.UI[TerrainNbX______].hasChanged()) isAllocated= false;
  if (D.UI[TerrainNbY______].hasChanged()) isAllocated= false;
  if (D.UI[DropletNbK______].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool TerrainErosion::CheckRefresh() {
  if (D.UI[TerrainNbCut____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void TerrainErosion::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  terrainNbX= std::max(2, D.UI[TerrainNbX______].I());
  terrainNbY= std::max(2, D.UI[TerrainNbY______].I());
  dropletNbK= std::max(1, D.UI[DropletNbK______].I());

  // Allocate data
  terrainPos= Field::Field2(terrainNbX, terrainNbY, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  terrainNor= Field::Field2(terrainNbX, terrainNbY, Vec::Vec3<float>(0.0f, 0.0f, 1.0f));
  terrainCol= Field::Field2(terrainNbX, terrainNbY, Vec::Vec3<float>(0.5f, 0.5f, 0.5f));
  terrainChg= Field::Field2(terrainNbX, terrainNbY, 0.0f);

  dropletPosOld= std::vector<Vec::Vec3<float>>(dropletNbK, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  dropletPosCur= std::vector<Vec::Vec3<float>>(dropletNbK, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  dropletVelCur= std::vector<Vec::Vec3<float>>(dropletNbK, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  dropletAccCur= std::vector<Vec::Vec3<float>>(dropletNbK, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  dropletForCur= std::vector<Vec::Vec3<float>>(dropletNbK, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  dropletColCur= std::vector<Vec::Vec3<float>>(dropletNbK, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  dropletMasCur= std::vector<float>(dropletNbK, 0.0f);
  dropletRadCur= std::vector<float>(dropletNbK, 0.0f);
  dropletSatCur= std::vector<float>(dropletNbK, 0.0f);
  dropletIsDead= std::vector<bool>(dropletNbK, true);
}


// Refresh the project
void TerrainErosion::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Get UI parameters
  terrainNbC= std::max(0, D.UI[TerrainNbCut____].I());

  // Reset random seed to always generate same terrain
  srand(0);

  // Precompute cut planes
  std::vector<Vec::Vec2<float>> cutPiv(terrainNbC);
  std::vector<Vec::Vec2<float>> cutVec(terrainNbC);
  for (int iter= 0; iter < terrainNbC; iter++) {
    cutPiv[iter].set(Random::Val(0.0f, 1.0f), Random::Val(0.0f, 1.0f));
    do {
      cutVec[iter].set(Random::Val(-1.0f, 1.0f), Random::Val(-1.0f, 1.0f));
    } while (cutVec[iter].normSquared() > 1.0f);
    cutVec[iter].normalize();
  }

  // Compute random terrain through iterative cutting
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      terrainPos.at(x, y)[0]= float(x) / float(terrainNbX - 1);
      terrainPos.at(x, y)[1]= float(y) / float(terrainNbY - 1);
      terrainPos.at(x, y)[2]= 0.0f;
      for (int iter= 0; iter < terrainNbC; iter++) {
        Vec::Vec2<float> pos(terrainPos.at(x, y)[0], terrainPos.at(x, y)[1]);
        if ((pos - cutPiv[iter]).dot(cutVec[iter]) < 0.0f)
          terrainPos.at(x, y)[2]+= 1.0f;
        else
          terrainPos.at(x, y)[2]-= 1.0f;
      }
    }
  }

  // Smooth the terrain
  Field::Field2<Vec::Vec3<float>> terrainPosOld= terrainPos;
  for (int iter= 0; iter < std::max(terrainNbX, terrainNbY) / 64; iter++) {
    terrainPosOld.swap(terrainPos);
    for (int x= 0; x < terrainNbX; x++) {
      for (int y= 0; y < terrainNbY; y++) {
        int count= 0;
        terrainPos.at(x, y)[2]= 0.0;
        for (int xOff= std::max(x - 1, 0); xOff <= std::min(x + 1, terrainNbX - 1); xOff++) {
          for (int yOff= std::max(y - 1, 0); yOff <= std::min(y + 1, terrainNbY - 1); yOff++) {
            terrainPos.at(x, y)[2]+= terrainPosOld.at(xOff, yOff)[2];
            count++;
          }
        }
        terrainPos.at(x, y)[2]/= float(count);
      }
    }
  }

  // Rescale terrain elevation
  float terrainMinTarg= 0.3f;
  float terrainMaxTarg= 0.7f;
  float terrainMinVal= terrainPos.at(0, 0)[2];
  float terrainMaxVal= terrainPos.at(0, 0)[2];
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      if (terrainMinVal > terrainPos.at(x, y)[2]) terrainMinVal= terrainPos.at(x, y)[2];
      if (terrainMaxVal < terrainPos.at(x, y)[2]) terrainMaxVal= terrainPos.at(x, y)[2];
    }
  }
  for (int x= 0; x < terrainNbX; x++)
    for (int y= 0; y < terrainNbY; y++)
      terrainPos.at(x, y)[2]= terrainMinTarg + (terrainMaxTarg - terrainMinTarg) * (terrainPos.at(x, y)[2] - terrainMinVal) / (terrainMaxVal - terrainMinVal);

  // Compute terrain mesh vertex normals
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      terrainNor.at(x, y).set(0.0f, 0.0f, 0.0f);
      if (x > 0 && y > 0)
        terrainNor.at(x, y)+= ((terrainPos.at(x - 1, y) - terrainPos.at(x, y)).cross(terrainPos.at(x, y - 1) - terrainPos.at(x, y))).normalized();
      if (x < terrainNbX - 1 && y > 0)
        terrainNor.at(x, y)+= ((terrainPos.at(x, y - 1) - terrainPos.at(x, y)).cross(terrainPos.at(x + 1, y) - terrainPos.at(x, y))).normalized();
      if (x < terrainNbX - 1 && y < terrainNbY - 1)
        terrainNor.at(x, y)+= ((terrainPos.at(x + 1, y) - terrainPos.at(x, y)).cross(terrainPos.at(x, y + 1) - terrainPos.at(x, y))).normalized();
      if (x > 0 && y < terrainNbY - 1)
        terrainNor.at(x, y)+= ((terrainPos.at(x, y + 1) - terrainPos.at(x, y)).cross(terrainPos.at(x - 1, y) - terrainPos.at(x, y))).normalized();
      terrainNor.at(x, y).normalize();
    }
  }
}


// Handle UI parameter change
void TerrainErosion::ParamChange() {
}


// Handle keypress
void TerrainErosion::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Handle mouse action
void TerrainErosion::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void TerrainErosion::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  float dt= D.UI[SimuTimestep____].F();
  float velocityDecay= std::min(std::max(D.UI[VelDecay________].F(), 0.0f), 1.0f);
  Vec::Vec3<float> gravity(0.0f, 0.0f, -0.5f);

  // Respawn dead droplets
  for (int k= 0; k < dropletNbK; k++) {
    if (dropletIsDead[k]) {
      dropletPosCur[k].set(Random::Val(0.0f, 1.0f), Random::Val(0.0f, 1.0f), Random::Val(0.7f, 1.0f));
      dropletPosOld[k]= dropletPosCur[k];
      dropletColCur[k].set(0.5f, 0.5f, 1.0f);
      dropletMasCur[k]= 1.0f;
      dropletSatCur[k]= 0.01f;
      dropletRadCur[k]= std::max(D.UI[DropletRad______].F(), 1.e-6f);
      dropletIsDead[k]= false;
    }
  }

  // Reset forces
  for (int k= 0; k < dropletNbK; k++)
    dropletForCur[k].set(0.0f, 0.0f, 0.0f);

  // Add gravity forces
  for (int k= 0; k < dropletNbK; k++)
    dropletForCur[k]+= gravity * dropletMasCur[k];

  // Reset terrain change
  for (int x= 0; x < terrainNbX; x++)
    for (int y= 0; y < terrainNbY; y++)
      terrainChg.at(x, y)= 0.0f;

  // Compute terrain erosion and sedimentation
  for (int k= 0; k < dropletNbK; k++) {
    int xRef= std::min(std::max(int(std::round(dropletPosCur[k][0] * float(terrainNbX - 1))), 0), terrainNbX - 1);
    int yRef= std::min(std::max(int(std::round(dropletPosCur[k][1] * float(terrainNbY - 1))), 0), terrainNbY - 1);
    int xRad= int(std::ceil(dropletRadCur[k] * 2.0f / (1.0f / float(terrainNbX))));
    int yRad= int(std::ceil(dropletRadCur[k] * 2.0f / (1.0f / float(terrainNbY))));
    for (int xOff= std::max(xRef - xRad, 0); xOff <= std::min(xRef + xRad, terrainNbX - 1); xOff++) {
      for (int yOff= std::max(yRef - yRad, 0); yOff <= std::min(yRef + yRad, terrainNbY - 1); yOff++) {
        float weight= std::max(dropletRadCur[k] * 2.0f - (dropletPosCur[k] - terrainPos.at(xOff, yOff)).norm(), 0.0f);
        terrainChg.at(xOff, yOff)-= weight * D.UI[ErosionCoeff____].F();
      }
    }
  }

  // Resolve droplet-terrain collisions
  for (int k= 0; k < dropletNbK; k++) {
    float xFloat= dropletPosCur[k][0] * float(terrainNbX - 1);
    float yFloat= dropletPosCur[k][1] * float(terrainNbY - 1);

    int x0= std::min(std::max(int(std::floor(xFloat)), 0), terrainNbX - 2);
    int y0= std::min(std::max(int(std::floor(yFloat)), 0), terrainNbY - 2);
    int x1= x0 + 1;
    int y1= y0 + 1;

    float xWeight1= xFloat - float(x0);
    float yWeight1= yFloat - float(y0);
    float xWeight0= 1.0 - xWeight1;
    float yWeight0= 1.0 - yWeight1;

    float interpoVal= 0.0;
    interpoVal+= terrainPos.at(x0, y0)[2] * (xWeight0 * yWeight0);
    interpoVal+= terrainPos.at(x0, y1)[2] * (xWeight0 * yWeight1);
    interpoVal+= terrainPos.at(x1, y0)[2] * (xWeight1 * yWeight0);
    interpoVal+= terrainPos.at(x1, y1)[2] * (xWeight1 * yWeight1);

    if (dropletPosCur[k][2] - dropletRadCur[k] < interpoVal) {
      Vec::Vec3<float> interpoNor(0.0f, 0.0f, 0.0f);
      interpoNor+= terrainNor.at(x0, y0) * (xWeight0 * yWeight0);
      interpoNor+= terrainNor.at(x0, y1) * (xWeight0 * yWeight1);
      interpoNor+= terrainNor.at(x1, y0) * (xWeight1 * yWeight0);
      interpoNor+= terrainNor.at(x1, y1) * (xWeight1 * yWeight1);
      dropletPosCur[k]+= (interpoVal + dropletRadCur[k] - dropletPosCur[k][2]) * interpoNor.normalized();
    }
  }

  // Resolve droplet-droplet collisions
  for (int k0= 0; k0 < dropletNbK; k0++) {
    for (int k1= k0 + 1; k1 < dropletNbK; k1++) {
      if ((dropletPosCur[k1] - dropletPosCur[k0]).normSquared() <= (dropletRadCur[k0] + dropletRadCur[k1]) * (dropletRadCur[k0] + dropletRadCur[k1])) {
        Vec::Vec3<float> val= (dropletPosCur[k1] - dropletPosCur[k0]).normalized() * 0.5f * ((dropletRadCur[k0] + dropletRadCur[k1]) - (dropletPosCur[k1] - dropletPosCur[k0]).norm());
        dropletPosCur[k0]-= val;
        dropletPosCur[k1]+= val;
      }
    }
  }

  // Deduce velocities
  for (int k= 0; k < dropletNbK; k++)
    dropletVelCur[k]= (dropletPosCur[k] - dropletPosOld[k]) / dt;

  // Apply explicit velocity damping
  for (int k= 0; k < dropletNbK; k++)
    dropletVelCur[k]= dropletVelCur[k] * std::pow(velocityDecay, dt);

  // Update positions
  dropletPosOld= dropletPosCur;
  for (int k= 0; k < dropletNbK; k++) {
    dropletAccCur[k]= dropletForCur[k] / dropletMasCur[k];
    dropletPosCur[k]= dropletPosCur[k] + dropletVelCur[k] * dt + dropletAccCur[k] * dt * dt;
  }

  // Kill particles out of bounds
  for (int k= 0; k < dropletNbK; k++) {
    if (dropletPosCur[k][0] < 0.0f || dropletPosCur[k][0] > 1.0f) dropletIsDead[k]= true;
    if (dropletPosCur[k][1] < 0.0f || dropletPosCur[k][1] > 1.0f) dropletIsDead[k]= true;
    if (dropletPosCur[k][2] < 0.0f || dropletPosCur[k][2] > 1.0f) dropletIsDead[k]= true;
  }

  // Apply terrain change
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      float minNeighbor= terrainPos.at(x, y)[2];
      float maxNeighbor= terrainPos.at(x, y)[2];
      for (int xOff= std::max(x - 1, 0); xOff <= std::min(x + 1, terrainNbX - 1); xOff++) {
        for (int yOff= std::max(y - 1, 0); yOff <= std::min(y + 1, terrainNbY - 1); yOff++) {
          if (minNeighbor > terrainPos.at(xOff, yOff)[2]) minNeighbor= terrainPos.at(xOff, yOff)[2];
          if (maxNeighbor < terrainPos.at(xOff, yOff)[2]) maxNeighbor= terrainPos.at(xOff, yOff)[2];
        }
      }
      terrainPos.at(x, y)[2]= std::min(std::max(terrainPos.at(x, y)[2] + terrainChg.at(x, y), minNeighbor), maxNeighbor);
    }
  }

  // Smooth the terrain
  Field::Field2<Vec::Vec3<float>> terrainPosOld= terrainPos;
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      int count= 0;
      terrainPos.at(x, y)[2]= 0.0;
      for (int xOff= std::max(x - 1, 0); xOff <= std::min(x + 1, terrainNbX - 1); xOff++) {
        for (int yOff= std::max(y - 1, 0); yOff <= std::min(y + 1, terrainNbY - 1); yOff++) {
          terrainPos.at(x, y)[2]+= terrainPosOld.at(xOff, yOff)[2];
          count++;
        }
      }
      terrainPos.at(x, y)[2]/= float(count);
      terrainPos.at(x, y)[2]= D.UI[SmoothResist____].F() * terrainPosOld.at(x, y)[2] + (1.0f - D.UI[SmoothResist____].F()) * terrainPos.at(x, y)[2];
    }
  }

  // Recompute terrain normals
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      terrainNor.at(x, y).set(0.0f, 0.0f, 0.0f);
      if (x > 0 && y > 0)
        terrainNor.at(x, y)+= ((terrainPos.at(x - 1, y) - terrainPos.at(x, y)).cross(terrainPos.at(x, y - 1) - terrainPos.at(x, y))).normalized();
      if (x < terrainNbX - 1 && y > 0)
        terrainNor.at(x, y)+= ((terrainPos.at(x, y - 1) - terrainPos.at(x, y)).cross(terrainPos.at(x + 1, y) - terrainPos.at(x, y))).normalized();
      if (x < terrainNbX - 1 && y < terrainNbY - 1)
        terrainNor.at(x, y)+= ((terrainPos.at(x + 1, y) - terrainPos.at(x, y)).cross(terrainPos.at(x, y + 1) - terrainPos.at(x, y))).normalized();
      if (x > 0 && y < terrainNbY - 1)
        terrainNor.at(x, y)+= ((terrainPos.at(x, y + 1) - terrainPos.at(x, y)).cross(terrainPos.at(x - 1, y) - terrainPos.at(x, y))).normalized();
      terrainNor.at(x, y).normalize();
    }
  }
}


// Draw the project
void TerrainErosion::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Set the terrain colors
  float terrainMinVal= terrainPos.at(0, 0)[2];
  float terrainMaxVal= terrainPos.at(0, 0)[2];
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      if (terrainMinVal > terrainPos.at(x, y)[2]) terrainMinVal= terrainPos.at(x, y)[2];
      if (terrainMaxVal < terrainPos.at(x, y)[2]) terrainMaxVal= terrainPos.at(x, y)[2];
    }
  }
  for (int x= 0; x < terrainNbX; x++) {
    for (int y= 0; y < terrainNbY; y++) {
      if (D.UI[ColorMode_______].I() == 0) {
        float val= (terrainPos.at(x, y)[2] - terrainMinVal) / (terrainMaxVal - terrainMinVal);
        Colormap::RatioToJetBrightSmooth(val, terrainCol.at(x, y)[0], terrainCol.at(x, y)[1], terrainCol.at(x, y)[2]);
      }
      else if (D.UI[ColorMode_______].I() == 1) {
        terrainCol.at(x, y)[0]= 0.5f + terrainNor.at(x, y)[0] / 2.0f;
        terrainCol.at(x, y)[1]= 0.5f + terrainNor.at(x, y)[1] / 2.0f;
        terrainCol.at(x, y)[2]= 0.5f + terrainNor.at(x, y)[2] / 2.0f;
      }
      else if (D.UI[ColorMode_______].I() == 2) {
        if (terrainNor.at(x, y).dot(Vec::Vec3<float>(0.0f, 0.0f, 1.0f)) < D.UI[CliffThresh_____].F()) {
          terrainCol.at(x, y)[0]= 0.7f;
          terrainCol.at(x, y)[1]= 0.6f;
          terrainCol.at(x, y)[2]= 0.3f;
        }
        else {
          terrainCol.at(x, y)[0]= 0.5f;
          terrainCol.at(x, y)[1]= 0.9f;
          terrainCol.at(x, y)[2]= 0.5f;
        }
      }
    }
  }

  // Draw the terrain
  if (D.displayMode[1]) {
    glEnable(GL_LIGHTING);
    glBegin(GL_QUADS);
    for (int x= 0; x < terrainNbX - 1; x++) {
      for (int y= 0; y < terrainNbY - 1; y++) {
        Vec::Vec3<float> flatNormal= (terrainNor.at(x, y) + terrainNor.at(x + 1, y) + terrainNor.at(x + 1, y + 1) + terrainNor.at(x, y + 1)).normalized();
        Vec::Vec3<float> flatColor= (terrainCol.at(x, y) + terrainCol.at(x + 1, y) + terrainCol.at(x + 1, y + 1) + terrainCol.at(x, y + 1)) / 4.0f;
        glColor3fv(flatColor.array());
        glNormal3fv(flatNormal.array());
        glVertex3fv(terrainPos.at(x, y).array());
        glVertex3fv(terrainPos.at(x + 1, y).array());
        glVertex3fv(terrainPos.at(x + 1, y + 1).array());
        glVertex3fv(terrainPos.at(x, y + 1).array());
      }
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }

  // Draw the terrain normals
  if (D.displayMode[4]) {
    glBegin(GL_LINES);
    for (int x= 0; x < terrainNbX; x++) {
      for (int y= 0; y < terrainNbY; y++) {
        float r, g, b;
        Colormap::RatioToJetSmooth(1.0f - terrainNor.at(x, y)[2] * terrainNor.at(x, y)[2], r, g, b);
        glColor3f(r, g, b);
        glVertex3fv(terrainPos.at(x, y).array());
        glVertex3fv((terrainPos.at(x, y) + 0.02f * terrainNor.at(x, y)).array());
      }
    }
    glEnd();
  }

  // Draw the droplets
  if (D.displayMode[5]) {
    for (int k= 0; k < dropletNbK; k++) {
      glPushMatrix();
      glTranslatef(dropletPosCur[k][0], dropletPosCur[k][1], dropletPosCur[k][2]);
      glScalef(dropletRadCur[k], dropletRadCur[k], dropletRadCur[k]);
      float r, g, b;
      Colormap::RatioToJetSmooth(dropletVelCur[k].norm(), r, g, b);
      glColor4f(r, g, b, 0.5f);
      glutSolidSphere(1.0, 32, 16);
      glPopMatrix();
    }
  }
}
