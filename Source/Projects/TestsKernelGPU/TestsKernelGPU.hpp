#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Test environment for running GPU kernel using OpenCL Wrapper
class TestsKernelGPU
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ArraySize_______,
    ______________00,
    NbParticles_____,
    InitVel_________,
    Timestep________,
    Epsilon_________,
    GravCoeff_______,
    ColorMode_______,
    ScaleColor______,
    ScaleShape______,
    ______________01,
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

  unsigned int N;
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  bool isOpenCLReady;

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
  void KeyPress(const unsigned char key);
  void MousePress(const unsigned char mouse);
  void Animate();
  void Draw();
};