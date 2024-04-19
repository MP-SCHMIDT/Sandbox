#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Improvised curved space-time rendering
// - 4D Eulerian representation of the world
// - Local curvature evaluated from mass distribution
// - Explicit integration of photons backtracing in the curved fields from the screen to the source
class SpaceTimeWorld
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    WorldNbT________,
    WorldNbX________,
    WorldNbY________,
    WorldNbZ________,
    ScreenNbH_______,
    ScreenNbV_______,
    ScreenNbS_______,
    CursorWorldT____,
    MassReach_______,
    TimePersist_____,
    FactorCurv______,
    FactorDoppl_____,
    VerboseLevel____,
  };

  int worldNbT;
  int worldNbX;
  int worldNbY;
  int worldNbZ;
  std::vector<std::vector<std::vector<std::vector<bool>>>> worldSolid;
  std::vector<std::vector<std::vector<std::vector<bool>>>> worldIsFix;
  std::vector<std::vector<std::vector<std::vector<float>>>> worldMasss;
  std::vector<std::vector<std::vector<std::vector<Vec::Vec3<float>>>>> worldColor;
  std::vector<std::vector<std::vector<std::vector<Vec::Vec4<float>>>>> worldFlows;

  int screenNbH;
  int screenNbV;
  int screenNbS;
  std::vector<std::vector<Vec::Vec3<float>>> screenColor;
  std::vector<std::vector<int>> screenCount;
  std::vector<std::vector<std::vector<Vec::Vec4<float>>>> photonPos;
  std::vector<std::vector<std::vector<Vec::Vec4<float>>>> photonVel;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  SpaceTimeWorld();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void Refresh();
  void KeyPress(const unsigned char key);
  void MousePress(const unsigned char mouse);
  void Animate();
  void Draw();
};
