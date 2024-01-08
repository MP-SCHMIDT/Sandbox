#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// String Art Generator
// - Algorithmic solver for the string art optimization problem
// - Reproduce the appearance of a target image by wrapping a string around a set of fixed pegs
class StringArtOptim
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ImageID_________,
    ImageSizeW______,
    ImageSizeH______,
    PegLayout_______,
    PegNumber_______,
    ColorsAdd_______,
    ColorsSub_______,
    StepCount_______,
    SingleLine______,
    BlendMode_______,
    CoeffColor______,
    VerboseLevel____,
  };

  // Problem dimensions
  int nW;
  int nH;

  // Fields for scenario setup
  std::vector<std::vector<Vec::Vec3<float>>> ImRef;
  std::vector<std::vector<Vec::Vec3<float>>> ImCur;
  std::vector<std::array<int, 2>> Pegs;
  std::vector<int> PegsCount;
  std::vector<std::vector<int>> Lines;
  std::vector<Vec::Vec3<float>> Colors;

  bool AddLineStep();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  StringArtOptim();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
