#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


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

  int mapNbX;
  int mapNbY;
  int mapNbIter;
  double mapDivThresh;
  double mapZoom;
  Vec::Vec2<double> mapFocus;
  Vec::Vec2<double> mapConst;

  std::vector<std::vector<Vec::Vec3<float>>> mapPos;
  std::vector<std::vector<Vec::Vec3<float>>> mapNor;
  std::vector<std::vector<Vec::Vec3<float>>> mapCol;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  FractalElevMap();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
