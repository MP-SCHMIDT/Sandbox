#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Vec.hpp"


// N-Body simulation
// - Barnes-Hut algorithm to account for gravitational interaction
//
// Reference
// https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation
// http://arborjs.org/docs/barnes-hut
class NBodyGravDynam
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    BodyCount_______,
    BodyRandSeed____,
    BodyInitLayout__,
    BodyInitVel_____,
    BodyRadius______,
    ______________00,
    DomainLock2D____,
    DomainTorusPos__,
    DomainTorusFor__,
    ______________01,
    TreeMaxDepth____,
    TreeInfiniteBox_,
    ______________02,
    SimuMode________,
    SimuTreeTol_____,
    SimuTimeStep____,
    SimuStepBatch___,
    SimuTotGravity__,
    SimuDrag________,
    SimuCollision___,
    SimuBodySort____,
    SimuMultithread_,
    ______________03,
    ColorMode_______,
    ColorFactor_____,
    ScaleFactor_____,
    SphereSimple____,
    ShowEmptyCells__,
    ______________04,
    TestParamNBS_00_,
    TestParamNBS_01_,
    TestParamNBS_02_,
    TestParamNBS_03_,
    TestParamNBS_04_,
    VerboseLevel____,
  };

  struct OctreeNode
  {
    float Size;                  // Size of the cell
    Vec::Vec3<float> Center;     // Center of the cell
    Vec::Vec3<float> AvgPos;     // Average Position of bodies in the cell
    Vec::Vec3<float> AvgVel;     // Average velocity of bodies in the cell
    unsigned int Count;          // Count of bodies in the cell
    unsigned int Child;          // Index of the first of the 8 contiguous children, order is X major and Z minor
    unsigned int Next;           // Index of the next node in the depth first search (ignoring the children)
  };

  unsigned int N;
  float simTime= 0.0f;

  double timerSort= 0.0;
  double timerTree= 0.0;
  double timerForces= 0.0;
  double timerSimu= 0.0;
  double timerDraw= 0.0;
  
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<Vec::Vec3<float>> For;
  std::vector<NBodyGravDynam::OctreeNode> Tree;
  // #define TESTING_DISPLAY_FORCES_VECTORS
  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  std::vector<Vec::Vec3<float>> ContribPos;
  std::vector<unsigned int> ContribCount;
  std::vector<unsigned int> ContribCell;
  #endif

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  NBodyGravDynam();

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

  private:
  // Display functions
  void DrawScene();

  // Simulation functions
  void StepSimulation();
  void ComputeForces(const std::vector<Vec::Vec3<float>>& iPos);
  void SetupGPU();
  void StepSimulationGPU();

  // Morton sort functions
  void MortonSortBodies();
  unsigned int MortonExpandBits(unsigned int v);
  unsigned int MortonGetIndex(float x, float y, float z);

  // Tree construction functions
  void BuildTree(const std::vector<Vec::Vec3<float>>& iPos);
  inline unsigned int PickSubCell(const Vec::Vec3<float> &iCenter,
                           const Vec::Vec3<float> &iPos,
                           const unsigned int iIdxFirstChild) {
    return iIdxFirstChild + (iPos[0] > iCenter[0]) * 4 + (iPos[1] > iCenter[1]) * 2 + (iPos[2] > iCenter[2]);
  }

  // Utility functions
  inline void UtilMake2D(Vec::Vec3<float> &ioPos, Vec::Vec3<float> &ioVel) {
    ioPos[0]= 0.49f;
    ioVel[0]= 0.0f;
  }
  inline void UtilMakeTorusPos(Vec::Vec3<float> &ioPos) {
    if (ioPos[0] < 0.0f || ioPos[0] > 1.0f) ioPos[0]-= std::floor(ioPos[0]);
    if (ioPos[1] < 0.0f || ioPos[1] > 1.0f) ioPos[1]-= std::floor(ioPos[1]);
    if (ioPos[2] < 0.0f || ioPos[2] > 1.0f) ioPos[2]-= std::floor(ioPos[2]);
  }
  inline void UtilMakeTorusFor(const Vec::Vec3<float> &iPos0, Vec::Vec3<float> &ioPos1) {
    if      ((ioPos1[0] - iPos0[0]) >  0.5f) ioPos1[0]-= 1.0f;
    else if ((ioPos1[0] - iPos0[0]) < -0.5f) ioPos1[0]+= 1.0f;
    if      ((ioPos1[1] - iPos0[1]) >  0.5f) ioPos1[1]-= 1.0f;
    else if ((ioPos1[1] - iPos0[1]) < -0.5f) ioPos1[1]+= 1.0f;
    if      ((ioPos1[2] - iPos0[2]) >  0.5f) ioPos1[2]-= 1.0f;
    else if ((ioPos1[2] - iPos0[2]) < -0.5f) ioPos1[2]+= 1.0f;
  }
};
