#pragma once

// Standard lib
#include <vector>

// Algo headers


// Generic 3D wave equation solver
// - ∂²u/∂t² = c² ∂²u/∂x²
// - Explicit integration
// - Centered finite difference discretization of operators
// https://vitalitylearning.medium.com/solving-the-1d-wave-equation-numerical-discretization-190a92c917bc
// https://en.wikipedia.org/wiki/Wave_equation
class WaveEquationFD
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ResolutionX_____,
    ResolutionY_____,
    ResolutionZ_____,
    VoxelSize_______,
    TimeStep________,
    MaxWaveSpeed____,
    ______________00,
    Multithread_____,
    ColorMode_______,
    ColorFactor_____,
    ______________01,
    TestParamWAV_0__,
    TestParamWAV_1__,
    TestParamWAV_2__,
    TestParamWAV_3__,
    TestParamWAV_4__,
    TestParamWAV_5__,
    TestParamWAV_6__,
    TestParamWAV_7__,
    TestParamWAV_8__,
    TestParamWAV_9__,
    ______________02,
    VerboseLevel____,
  };

  int nX;
  int nY;
  int nZ;
  double simTime;

  std::vector<std::vector<std::vector<double>>> UNew;   // Field value at new time
  std::vector<std::vector<std::vector<double>>> UCur;   // Field value at current time
  std::vector<std::vector<std::vector<double>>> UOld;   // Field value at previous time
  std::vector<std::vector<std::vector<double>>> Speed;  // Wave propagation speed


  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  WaveEquationFD();

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
