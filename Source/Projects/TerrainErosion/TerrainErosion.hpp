#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Field.hpp"
#include "Type/Vec.hpp"


// Terrain generation and erosion simulation
// - Representation as height map
// - Initial terrain created with iterative random cut planes
// - Particles dropped and collide with the terrain using Position Based Dynamics scheme
// - Erosion and sedimentation handled by heuristic rules
//
// Reference
// https://www.youtube.com/watch?v=eaXk97ujbPQ
class TerrainErosion
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    TerrainNbX______,
    TerrainNbY______,
    TerrainNbCut____,
    DropletNbK______,
    DropletRad______,
    SimuTimestep____,
    VelDecay________,
    ErosionCoeff____,
    SmoothResist____,
    CliffThresh_____,
    ColorMode_______,
    VerboseLevel____,
  };

  int terrainNbX;
  int terrainNbY;
  int terrainNbC;
  Field::Field2<Vec::Vec3<float>> terrainPos;
  Field::Field2<Vec::Vec3<float>> terrainNor;
  Field::Field2<Vec::Vec3<float>> terrainCol;
  Field::Field2<float> terrainChg;

  int dropletNbK;
  std::vector<Vec::Vec3<float>> dropletPosOld;
  std::vector<Vec::Vec3<float>> dropletPosCur;
  std::vector<Vec::Vec3<float>> dropletVelCur;
  std::vector<Vec::Vec3<float>> dropletAccCur;
  std::vector<Vec::Vec3<float>> dropletForCur;
  std::vector<Vec::Vec3<float>> dropletColCur;
  std::vector<float> dropletMasCur;
  std::vector<float> dropletRadCur;
  std::vector<float> dropletSatCur;
  std::vector<bool> dropletIsDead;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  TerrainErosion();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void Refresh();
  void ParamChange();
  void KeyPress();
  void MousePress();
  void Animate();
  void Draw();
};
