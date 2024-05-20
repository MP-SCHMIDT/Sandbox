#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers


// Generic 3D wave equation solver
// - ∂²u/∂t² = c² ∂²u/∂x²
// - Explicit integration
// - Centered finite difference discretization of operators
// - Reflection BC with Dirichlet
// - Absorption BC with Perfectly Matched Layer
// - Spatially varying wave speed
class WavePropagSimu
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ScenarioPreset__,
    ScenarioFileID__,
    ResolutionX_____,
    ResolutionY_____,
    ResolutionZ_____,
    VoxelSize_______,
    ______________00,
    TimeStep________,
    Parallelize_____,
    NbSubsteps______,
    MaxWaveSpeed____,
    MaxAmplitude____,
    SourceRadius____,
    SourceFrequency_,
    ______________01,
    SliceDim________,
    SliceRelPosX____,
    SliceRelPosY____,
    SliceRelPosZ____,
    ColorMode_______,
    ColorFactor_____,
    ScaleFactor_____,
    AlphaFactor_____,
    ______________02,
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
    ______________03,
    VerboseLevel____,
  };

  int nX;
  int nY;
  int nZ;
  double simTime;
  std::array<double, 3> sourcePos;

  std::vector<std::vector<std::vector<double>>> UNew;   // Field value at new time
  std::vector<std::vector<std::vector<double>>> UCur;   // Field value at current time
  std::vector<std::vector<std::vector<double>>> UOld;   // Field value at previous time
  std::vector<std::vector<std::vector<double>>> Speed;  // Wave propagation speed

  void StepSimulation();
  void AddSource();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  WavePropagSimu();

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
