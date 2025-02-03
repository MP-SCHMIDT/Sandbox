#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Vec.hpp"


// N-Body simulation
// - Barnes-Hut algorithm to account for gravitational interaction
//
// Reference
// todo
class NBodyGravDynam
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    BodyCount_______,
    BodyStartLayout_,
    BodyInitVel_____,
    BodyRadius______,
    ______________00,
    Domain2D________,
    DomainTaurusPos_,
    DomainTaurusFor_,
    ______________01,
    TreeMaxDepth____,
    TreeTolRatio____,
    TreeShowMin_____,
    TreeShowMax_____,
    TreeShowEmpty___,
    ______________02,
    SimuStepMode____,
    SimuForceMode___,
    SimuTimeStep____,
    SimuTotGravity__,
    SimuVelDecay____,
    ______________03,
    ColorMode_______,
    ColorFactor_____,
    ______________04,
    TestParamNBS_00_,
    TestParamNBS_01_,
    TestParamNBS_02_,
    TestParamNBS_03_,
    TestParamNBS_04_,
    TestParamNBS_05_,
    TestParamNBS_06_,
    TestParamNBS_07_,
    TestParamNBS_08_,
    TestParamNBS_09_,
    VerboseLevel____,
  };

  struct OctreeNode
  {
    Vec::Vec3<float> Center;     // Center of the cell
    float Size;                  // Size of the cell
    unsigned char Depth;         // Depth of the current cell in the tree
    Vec::Vec3<float> CenterMass; // Position of the center of mass in the cell
    unsigned int Count;          // Count of bodies in the cell
    unsigned int Child;          // Index of the first of the 8 contiguoius children, order is X major and Z minor
  };

  unsigned int N;
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<Vec::Vec3<float>> Acc;
  std::vector<Vec::Vec3<float>> For;
  std::vector<NBodyGravDynam::OctreeNode> Tree;
  #define TESTING_DISPLAY_FORCES_VECTORS
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
  void AddTreeNode(const float iCenterX,
                   const float iCenterY,
                   const float iCenterZ,
                   const float iHalfSize,
                   const unsigned char iDepth);
  void BuildTree();
  void ComputeForces();
  void ApplyPosBC();
};
