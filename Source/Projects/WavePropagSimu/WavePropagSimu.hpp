#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Type/Field.hpp"


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

  Field::Field3<double> UNew;                         // Field value at new time
  Field::Field3<double> UCur;                         // Field value at current time
  Field::Field3<double> UOld;                         // Field value at previous time
  Field::Field3<double> Speed;                        // Wave propagation speed
  Field::Field3<std::array<char, 3>> BorderNeighbor;  // Index offsets to neighbor if border cell

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
  void ParamChange();
  void KeyPress();
  void MousePress();
  void Animate();
  void Draw();
};
