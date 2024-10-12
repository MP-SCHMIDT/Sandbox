#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Field.hpp"
#include "Type/Vec.hpp"


// Generate 2.5D surface based on the Julia Set fractal
// - Representation with height map
// - Interactive depth, coefficients and zoom factor as UI parameters
//
// https://en.wikipedia.org/wiki/Julia_set
class FractalElevMap
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ResolutionX_____,
    ResolutionY_____,
    ZoomLevel_______,
    NumItera________,
    ShiftX__________,
    ShiftY__________,
    CoeffA__________,
    CoeffB__________,
    HeightLvls______,
    VerboseLevel____,
  };

  int nX;
  int nY;
  int mapNbIter;
  double mapDivThresh;
  double mapZoom;
  Vec::Vec2<double> mapFocus;
  Vec::Vec2<double> mapConst;

  Field::Field2<Vec::Vec3<float>> mapPos;
  Field::Field2<Vec::Vec3<float>> mapNor;
  Field::Field2<Vec::Vec3<float>> mapCol;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  FractalElevMap();

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
