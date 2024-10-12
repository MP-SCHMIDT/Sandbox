#pragma once


// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Type/Field.hpp"


class ImageExtruMesh
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ResolutionX_____,
    ResolutionY_____,
    ResolutionZ_____,
    SizeX___________,
    SizeY___________,
    SizeZ___________,
    CenterX_________,
    CenterY_________,
    CenterZ_________,
    HeightMapMode___,
    HeightMapSlope__,
    ValOffset_______,
    BaseRelHeight___,
    GeomSnapMode____,
    MidRelHeight____,
    SphereShift_____,
    OffsetScaling___,
    SmoothIter______,
    Isovalue________,
    VerboseLevel____,
  };

  Field::Field3<double> ScalarField;

  std::vector<std::array<double, 3>> Verts;
  std::vector<std::array<double, 3>> VertsCol;
  std::vector<std::array<int, 3>> Tris;
  std::vector<std::array<int, 4>> Quads;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  ImageExtruMesh();

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
