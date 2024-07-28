#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Field.hpp"


// Skeleton Folder
// - Empty template project
class AlgoTestEnviro
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    TestParamALG_00_,
    TestParamALG_01_,
    TestParamALG_02_,
    TestParamALG_03_,
    TestParamALG_04_,
    TestParamALG_05_,
    TestParamALG_06_,
    TestParamALG_07_,
    TestParamALG_08_,
    TestParamALG_09_,
    TestParamALG_10_,
    TestParamALG_11_,
    TestParamALG_12_,
    TestParamALG_13_,
    TestParamALG_14_,
    TestParamALG_15_,
    TestParamALG_16_,
    TestParamALG_17_,
    TestParamALG_18_,
    TestParamALG_19_,
    Isocut__________,
    ColorFactor_____,
    VerboseLevel____,
  };

  Field::Field3<double> ScalarField;
  std::vector<std::array<double, 3>> Verts;
  std::vector<std::array<int, 2>> Bars;
  std::vector<std::array<int, 3>> Tris;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  AlgoTestEnviro();

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
