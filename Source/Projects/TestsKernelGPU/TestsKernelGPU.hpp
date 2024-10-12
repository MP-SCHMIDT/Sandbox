#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Vec.hpp"


// Test environment for running GPU kernel using OpenCL Wrapper
// https://github.com/ProjectPhysX/OpenCL-Wrapper
class TestsKernelGPU
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ArraySize_______,
    ______________00,
    ReducSumSize____,
    ______________01,
    NbParticles_____,
    InitVel_________,
    Timestep________,
    Epsilon_________,
    GravCoeff_______,
    ColorMode_______,
    ScaleColor______,
    ScaleShape______,
    ______________02,
    TestParamGPU_00_,
    TestParamGPU_01_,
    TestParamGPU_02_,
    TestParamGPU_03_,
    TestParamGPU_04_,
    TestParamGPU_05_,
    TestParamGPU_06_,
    TestParamGPU_07_,
    TestParamGPU_08_,
    TestParamGPU_09_,
    VerboseLevel____,
  };

  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;

  void RunVecAddGPU();
  void RunReducSumGPU();
  void StepNBodySimGPU();
  void StepNBodySimCPU();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  TestsKernelGPU();

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
