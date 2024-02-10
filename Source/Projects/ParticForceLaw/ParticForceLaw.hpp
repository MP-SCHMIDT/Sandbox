#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Particle-based physics simulation
// - Memoryless point particles defined only by position and velocity
// - Different from DEM particles because no history, no rotation, etc.
// - Explicit time integration with forward euler or velocity verlet
// - Isotropic force law defining force vs distance
//
// Reference article from MIT CBA group
// Mesoscale material modeling with memoryless isotropic point particles
// https://www.sciencedirect.com/science/article/pii/S1877750323002582
class ParticForceLaw
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    DomainX_________,
    DomainY_________,
    DomainZ_________,
    LatticePitch____,
    LatticePattern__,
    ScenarioPreset__,
    ScenarioForce___,
    StepsPerDraw____,
    TimeStep________,
    IntegType_______,
    ParticleMass____,
    VelocityDamping_,
    ForceLawPreset__,
    ForceLawNormali_,
    ForceLawScale___,
    ForceLawCoeff00_,
    ForceLawCoeff01_,
    ForceLawCoeff02_,
    ForceLawCoeff03_,
    ForceLawCoeff04_,
    ForceLawCoeff05_,
    ForceLawCoeff06_,
    ForceLawCoeff07_,
    ForceLawCoeff08_,
    ForceLawCoeff09_,
    ForceLawCoeff10_,
    ForceLawCoeff11_,
    ForceLawCoeff12_,
    ForceLawCoeff13_,
    ForceLawCoeff14_,
    ColorMode_______,
    ColorFactor_____,
    VisuScale_______,
    VerboseLevel____,
  };

  // Force law
  std::vector<float> ForceLaw;

  // Particles
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<Vec::Vec3<float>> Acc;
  std::vector<Vec::Vec3<float>> For;
  std::vector<Vec::Vec3<float>> Col;
  std::vector<Vec::Vec3<float>> Ext;
  std::vector<int> Fix;

  void ComputeForces();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  ParticForceLaw();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
