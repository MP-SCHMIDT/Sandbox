#include "SpaceTimeWorld.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/Bresenham.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


class Shape
{
  public:
  int type;
  Vec::Vec3<float> posBeg;
  Vec::Vec3<float> posEnd;
  Vec::Vec3<float> col;
  float mass;
  float rad0;
  float rad1;

  Shape(
      int const iType,
      Vec::Vec3<float> const iPosBeg,
      Vec::Vec3<float> const iPosEnd,
      Vec::Vec3<float> const iCol,
      float const iMass,
      float const iRad0,
      float const iRad1) {
    type= iType;
    posBeg= iPosBeg;
    posEnd= iPosEnd;
    col= iCol;
    mass= iMass;
    rad0= iRad0;
    rad1= iRad1;
  }

  float ImplicitEval(Vec::Vec3<float> const iPosCell, float const iTimeRatio) {
    Vec::Vec3<float> posObj= posBeg + (posEnd - posBeg) * iTimeRatio;

    float val= 1.0f;
    if (type == 0)
      val= (posObj - iPosCell).norm() - rad0;
    if (type == 1)
      val= std::pow(std::sqrt(std::pow((posObj - iPosCell)[1], 2.0) + std::pow((posObj - iPosCell)[2], 2.0)) - rad0, 2.0) + std::pow((posObj - iPosCell)[0], 2.0) - std::pow(rad1, 2.0);

    return val;
  }
};


void SpaceTimeWorld::InitVoxelWorld() {
  // Get UI parameters
  worldNbT= std::max(D.UI[WorldNbT________].I(), 1);
  worldNbX= std::max(D.UI[WorldNbX________].I(), 1);
  worldNbY= std::max(D.UI[WorldNbY________].I(), 1);
  worldNbZ= std::max(D.UI[WorldNbZ________].I(), 1);

  // Allocate data
  worldSolid= Field::AllocNested4(worldNbT, worldNbX, worldNbY, worldNbZ, false);
  worldIsFix= Field::AllocNested4(worldNbT, worldNbX, worldNbY, worldNbZ, false);
  worldMasss= Field::AllocNested4(worldNbT, worldNbX, worldNbY, worldNbZ, 0.0f);
  worldColor= Field::AllocNested4(worldNbT, worldNbX, worldNbY, worldNbZ, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  worldFlows= Field::AllocNested4(worldNbT, worldNbX, worldNbY, worldNbZ, Vec::Vec4<float>(0.0f, 0.0f, 0.0f, 0.0f));

  // Load the BMP image for the background
  static std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (imageRGBA.empty()) {
    if (D.UI[InputFileID_____].I() == 0) FileInput::LoadImageBMPFile("./FileInput/Images/DeepField.bmp", imageRGBA, false);
    else if (D.UI[InputFileID_____].I() == 1) FileInput::LoadImageBMPFile("./FileInput/Images/PillarsCreation.bmp", imageRGBA, false);
    else FileInput::LoadImageBMPFile("./FileInput/Images/AlbertArt.bmp", imageRGBA, false);
  }

  // Add the background
  for (int t= 0; t < worldNbT; t++) {
    for (int x= 0; x < worldNbX; x++) {
      for (int y= 0; y < worldNbY; y++) {
        for (int z= 0; z < worldNbZ; z++) {
          // Add boundary conditions to spatial domain border
          if (x % worldNbX == 0 || y % worldNbY == 0 || z % worldNbZ == 0) {
            worldIsFix[t][x][y][z]= true;
            worldMasss[t][x][y][z]= 0.0f;
          }

          // Add grid background layer
          if (x == 0) {
            worldSolid[t][x][y][z]= true;
            worldColor[t][x][y][z].set(0.6f, 0.6f, 0.6f);
            if ((y + 1) % 8 <= 1 || (z + 1) % 8 <= 1)
              worldColor[t][x][y][z].set(0.4f, 0.4f, 0.4f);
          }
          // if (y == 0 || y == worldNbY - 1) {
          //   worldSolid[t][x][y][z]= true;
          //   worldColor[t][x][y][z].set(0.6f, 0.6f, 0.6f);
          //   if ((x + 1) % 10 <= 1 || (z + 1) % 10 <= 1)
          //     worldColor[t][x][y][z].set(0.4f, 0.4f, 0.4f);
          // }
          // if (z == 0 || z == worldNbZ - 1) {
          //   worldSolid[t][x][y][z]= true;
          //   worldColor[t][x][y][z].set(0.6f, 0.6f, 0.6f);
          //   if ((x + 1) % 10 <= 1 || (y + 1) % 10 <= 1)
          //     worldColor[t][x][y][z].set(0.4f, 0.4f, 0.4f);
          // }

          // Add image background layer
          if (x == 0) {
            int imgY= y * int(imageRGBA.size()) / worldNbY;
            int imgZ= z * int(imageRGBA[0].size()) / worldNbZ;
            worldSolid[t][x][y][z]= true;
            worldColor[t][x][y][z].set(imageRGBA[imgY][imgZ][0], imageRGBA[imgY][imgZ][1], imageRGBA[imgY][imgZ][2]);
          }
        }
      }
    }
  }

  // Create list of shapes to add
  std::vector<Shape> shapes;
  shapes.push_back(Shape(0, Vec::Vec3<float>(0.6f, -0.2f, -0.2f), Vec::Vec3<float>(0.6f, +1.2f, 1.2f), Vec::Vec3<float>(0.2f, 0.6f, 0.2f), +10.0f, 0.05f, 0.00f));  // 2 crossing balls
  // shapes.push_back(Shape(0, Vec::Vec3<float>(0.3f, +1.2f, -0.2f), Vec::Vec3<float>(0.3f, -0.2f, 1.2f), Vec::Vec3<float>(0.6f, 0.2f, 0.2f), -10.0f, 0.05f, 0.00f));  // 2 crossing balls
  // shapes.push_back(Shape(0, Vec::Vec3<float>(0.2f, +0.5f, +0.5f), Vec::Vec3<float>(0.8f, +0.5f, 0.5f), Vec::Vec3<float>(0.6f, 0.2f, 0.2f), +20.0f, 0.03f, 0.00f));  // 1 small approaching ball
  // shapes.push_back(Shape(1, Vec::Vec3<float>(0.5f, -0.5f, +0.5f), Vec::Vec3<float>(0.5f, +1.5f, 0.5f), Vec::Vec3<float>(0.3f, 0.3f, 0.7f), +10.0f, 0.15f, 0.03f));  // 1 moving donut

  // Add the shapes
  for (int t= 0; t < worldNbT; t++) {
    float posT= 0.5f;
    if (worldNbT > 1) posT= float(t) / float(worldNbT - 1);
    for (int x= 0; x < worldNbX; x++) {
      for (int y= 0; y < worldNbY; y++) {
        for (int z= 0; z < worldNbZ; z++) {
          Vec::Vec3<float> posCell((float(x) + 0.5f) / float(worldNbX), (float(y) + 0.5f) / float(worldNbY), (float(z) + 0.5f) / float(worldNbZ));
          for (Shape shape : shapes) {
            if (shape.ImplicitEval(posCell, posT) < 0.0f) {
              worldSolid[t][x][y][z]= true;
              worldIsFix[t][x][y][z]= true;
              worldMasss[t][x][y][z]= shape.mass;
              worldColor[t][x][y][z]= shape.col;
              // worldColor[t][x][y][z]= (1.0f - 0.6f * (posT-posCell)[1] / shape.rad0) * shape.col;
              // worldColor[t][x][y][z]= ((x + y + z) % 2 == 0) ? shape.col : 0.8f * shape.col;
            }
          }
        }
      }
    }
  }
}


void SpaceTimeWorld::ComputeWorldFlow() {
  // Precompute a mask for the world flow
  int maskSize= D.UI[MassReach_______].I();
  std::vector<std::vector<std::vector<std::vector<Vec::Vec4<float>>>>> maskVec= Field::AllocNested4(2 * maskSize + 1, 2 * maskSize + 1, 2 * maskSize + 1, 2 * maskSize + 1, Vec::Vec4<float>(0.0f, 0.0f, 0.0f, 0.0f));
  for (int t= 0; t < maskSize * 2 + 1; t++) {
    for (int x= 0; x < maskSize * 2 + 1; x++) {
      for (int y= 0; y < maskSize * 2 + 1; y++) {
        for (int z= 0; z < maskSize * 2 + 1; z++) {
          if (t == maskSize && x == maskSize && y == maskSize && z == maskSize) continue;
          Vec::Vec4<float> vec(float(maskSize - t), float(maskSize - x), float(maskSize - y), float(maskSize - z));
          maskVec[t][x][y][z]= vec.normalized() / vec.normSquared();
        }
      }
    }
  }

  // Compute the world flow
#pragma omp parallel for
  for (int t= 0; t < worldNbT; t++)
    for (int x= 0; x < worldNbX; x++)
      for (int y= 0; y < worldNbY; y++)
        for (int z= 0; z < worldNbZ; z++)
          worldFlows[t][x][y][z]= Vec::Vec4<float>(0.0f, 0.0f, 0.0f, 0.0f);

#pragma omp parallel for
  for (int t= 0; t < worldNbT; t++) {
    for (int x= 0; x < worldNbX; x++) {
      for (int y= 0; y < worldNbY; y++) {
        for (int z= 0; z < worldNbZ; z++) {
          if (worldMasss[t][x][y][z] == 0.0) continue;
          // for (int tOff= t; tOff <= t; tOff++)
          for (int tOff= t; tOff <= std::min(t + maskSize, worldNbT - 1); tOff++)
            // for (int tOff= std::max(t - maskSize, 0); tOff <= std::min(t + maskSize, worldNbT - 1); tOff++)
            for (int xOff= std::max(x - maskSize, 0); xOff <= std::min(x + maskSize, worldNbX - 1); xOff++)
              for (int yOff= std::max(y - maskSize, 0); yOff <= std::min(y + maskSize, worldNbY - 1); yOff++)
                for (int zOff= std::max(z - maskSize, 0); zOff <= std::min(z + maskSize, worldNbZ - 1); zOff++)
                  worldFlows[tOff][xOff][yOff][zOff]+= worldMasss[t][x][y][z] * maskVec[maskSize + tOff - t][maskSize + xOff - x][maskSize + yOff - y][maskSize + zOff - z];
        }
      }
    }
  }

  // // Add persistance between timesteps
  // for (int t= 1; t < worldNbT; t++)
  //   for (int x= 0; x < worldNbX; x++)
  //     for (int y= 0; y < worldNbY; y++)
  //       for (int z= 0; z < worldNbZ; z++)
  //         worldFlows[t][x][y][z]= worldFlows[t][x][y][z] + D.UI[TimePersist_].Get() * worldFlows[t - 1][x][y][z];
}


void SpaceTimeWorld::ComputeScreen() {
  // Get UI parameters
  screenNbH= std::max(D.UI[ScreenNbH_______].I(), 1);
  screenNbV= std::max(D.UI[ScreenNbV_______].I(), 1);
  screenNbS= std::max(D.UI[ScreenNbS_______].I(), 1);

  // Allocate data
  screenColor= Field::AllocNested2(screenNbH, screenNbV, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  screenCount= Field::AllocNested2(screenNbH, screenNbV, 1);
  photonPos= Field::AllocNested3(screenNbH, screenNbV, screenNbS, Vec::Vec4<float>(0.0f, 0.0f, 0.0f, 0.0f));
  photonVel= Field::AllocNested3(screenNbH, screenNbV, screenNbS, Vec::Vec4<float>(0.0f, 0.0f, 0.0f, 0.0f));

  // Initialize the photon fields
  for (int h= 0; h < screenNbH; h++) {
    for (int v= 0; v < screenNbV; v++) {
      screenColor[h][v]= Vec::Vec3<float>(0.0f, 0.0f, 0.0f);
      screenCount[h][v]= 1;
      photonPos[h][v][0]= Vec::Vec4<float>(((float)std::min(std::max(D.UI[CursorWorldT____].I(), 0), worldNbT - 1) + 0.5f) / float(worldNbT),
                                           1.0f - 1.0f / float(worldNbX),
                                           (0.5f + float(h)) / float(screenNbH),
                                           (0.5f + float(v)) / float(screenNbV));
      photonVel[h][v][0]= Vec::Vec4<float>(0.0f, -2.0f, 0.0f, 0.0f);
    }
  }

  // Compute the photon paths to render the scene on the screen
#pragma omp parallel for
  for (int h= 0; h < screenNbH; h++) {
    for (int v= 0; v < screenNbV; v++) {
      for (int s= 0; s < screenNbS - 1; s++) {
        int idxT= int(std::floor(photonPos[h][v][s][0] * float(worldNbT)));
        int idxX= int(std::floor(photonPos[h][v][s][1] * float(worldNbX)));
        int idxY= int(std::floor(photonPos[h][v][s][2] * float(worldNbY)));
        int idxZ= int(std::floor(photonPos[h][v][s][3] * float(worldNbZ)));
        if (idxT < 0 || idxT >= worldNbT || idxX < 0 || idxX >= worldNbX || idxY < 0 || idxY >= worldNbY || idxZ < 0 || idxZ >= worldNbZ) break;

        // if (idxT < 0 || idxT >= worldNbT || idxX < 0 || idxX >= worldNbX || idxY < 0 || idxY >= worldNbY || idxZ < 0 || idxZ >= worldNbZ) {
        //   float velDif= D.UI[FactorDoppl_].Get() * (photonVel[h][v][0].norm() - photonVel[h][v][s].norm());
        //   screenColor[h][v]= Vec::Vec3<float>(0.1, 0.1, 0.1) * (1.0 + velDif);
        //   break;
        // }

        photonPos[h][v][s + 1]= photonPos[h][v][s] + photonVel[h][v][s] / float(screenNbS);
        Vec::Vec4<float> flow(
            worldFlows[idxT][idxX][idxY][idxZ][0] * D.UI[TimePersist_____].F(),
            worldFlows[idxT][idxX][idxY][idxZ][1] * D.UI[FactorCurv______].F(),
            worldFlows[idxT][idxX][idxY][idxZ][2] * D.UI[FactorCurv______].F(),
            worldFlows[idxT][idxX][idxY][idxZ][3] * D.UI[FactorCurv______].F());
        photonVel[h][v][s + 1]= photonVel[h][v][s] + D.UI[FactorCurv______].F() * flow / float(screenNbS);
        screenCount[h][v]++;

        int endX= std::min(std::max(int(std::floor(photonPos[h][v][s + 1][1] * worldNbX)), 0), worldNbX - 1);
        int endY= std::min(std::max(int(std::floor(photonPos[h][v][s + 1][2] * worldNbY)), 0), worldNbY - 1);
        int endZ= std::min(std::max(int(std::floor(photonPos[h][v][s + 1][3] * worldNbZ)), 0), worldNbZ - 1);
        bool foundCollision= false;
        std::vector<std::array<int, 3>> listVox= Bresenham::Line3D(idxX, idxY, idxZ, endX, endY, endZ);
        for (std::array<int, 3> voxIdx : listVox) {
          if (worldSolid[idxT][voxIdx[0]][voxIdx[1]][voxIdx[2]]) {
            // float velDif= D.UI[FactorDoppl_].Get() * (photonVel[h][v][0].norm() - photonVel[h][v][s].norm());
            // screenColor[h][v]= worldColor[idxT][voxIdx[0]][voxIdx[1]][voxIdx[2]] * (1.0 + velDif);
            float velDif= D.UI[FactorDoppl_____].F() * (photonPos[h][v][0][0] - photonPos[h][v][s][0]);
            screenColor[h][v]= worldColor[idxT][voxIdx[0]][voxIdx[1]][voxIdx[2]] * (1.0f + velDif);
            foundCollision= true;
            break;
          }
        }
        // if (worldSolid[idxT][endX][endY][endZ]) {
        //   // float velDif= D.UI[FactorDoppl_].Get() * (photonVel[h][v][0].norm() - photonVel[h][v][s].norm());
        //   // screenColor[h][v]= worldColor[idxT][voxIdx[0]][voxIdx[1]][voxIdx[2]] * (1.0f + velDif);
        //   float velDif= D.UI[FactorDoppl_].Get() * (photonPos[h][v][0][0] - photonPos[h][v][s][0]);
        //   screenColor[h][v]= worldColor[idxT][endX][endY][endZ] * (1.0f + velDif);
        //   foundCollision= true;
        //   break;
        // }
        if (foundCollision) break;
      }
    }
  }
}
