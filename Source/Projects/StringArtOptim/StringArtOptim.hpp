#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Field.hpp"
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
    ColorsNormalize_,
    StepCount_______,
    SingleLine______,
    BlendMode_______,
    WeightOver______,
    WeightUnder_____,
    CoeffColor______,
    MinRelChange____,
    VerboseLevel____,
  };

  // Problem dimensions
  int nW;
  int nH;

  // Fields for scenario setup
  Field::Field2<Vec::Vec3<float>> ImRef;  // RGB reference image
  Field::Field2<Vec::Vec3<float>> ImCur;  // RGB current image
  std::vector<std::array<int, 2>> Pegs;   // 2D peg positions in pixel space
  std::vector<int> PegsUseCount;          // Counter for the number of uses for each peg
  std::vector<std::vector<int>> Stings;   // Path in peg ID for the string of each color
  std::vector<Vec::Vec3<float>> Colors;   // Colors of the strings

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
  void Refresh();
  void ParamChange();
  void KeyPress();
  void MousePress();
  void Animate();
  void Draw();
};
