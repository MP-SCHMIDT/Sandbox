#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers


// Skeleton Folder
// - Empty template project
class AlgoTestEnviro
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    TestParam00_____,
    TestParam01_____,
    TestParam02_____,
    TestParam03_____,
    TestParam04_____,
    TestParam05_____,
    TestParam06_____,
    TestParam07_____,
    TestParam08_____,
    TestParam09_____,
    TestParam10_____,
    TestParam11_____,
    TestParam12_____,
    TestParam13_____,
    TestParam14_____,
    TestParam15_____,
    TestParam16_____,
    TestParam17_____,
    TestParam18_____,
    TestParam19_____,
    Isocut__________,
    ColorFactor_____,
    VerboseLevel____,
  };

  std::vector<std::vector<std::vector<double>>> ScalarField;
  std::vector<std::vector<std::vector<std::array<double, 3>>>> VectorField;
  std::vector<std::vector<std::vector<std::array<double, 9>>>> TensorField;

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
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
