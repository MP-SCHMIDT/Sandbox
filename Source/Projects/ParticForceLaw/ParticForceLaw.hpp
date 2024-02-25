#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Particle-based physics simulation
// - Memoryless Isotropic Point Particles (MIPP) defined only by position and velocity
// - Different from DEM simulations because no history, no rotation, etc.
// - Base particle cloud as FCC lattice, Poisson sphere sampling, etc.
// - Explicit time integration with forward Euler or velocity Verlet
// - Isotropic force law defining force vs distance
// - Fast neighborhood search with spatial partition into buckets of fixed size
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
    ScenarioPreset__,
    Scenario2DID____,
    Scenario2DThick_,
    LatticePitch____,
    LatticePattern__,
    ______________00,
    ForceLawPreset__,
    ForceLawNormali_,
    ForceLawScale___,
    ForceLawA_______,
    ForceLaw08______,
    ForceLaw09______,
    ForceLawB_______,
    ForceLaw11______,
    ForceLaw12______,
    ForceLaw13______,
    ForceLawC_______,
    ForceLaw15______,
    ForceLaw20______,
    ForceLaw25______,
    ForceLaw30______,
    ______________01,
    BCVelX__________,
    BCVelY__________,
    BCVelZ__________,
    BCForX__________,
    BCForY__________,
    BCForZ__________,
    StepsPerDraw____,
    TimeStep________,
    MaterialDensity_,
    RadialDamping___,
    VelocityDamping_,
    BucketsCount____,
    BucketsCapacity_,
    IntegType_______,
    UseForceControl_,
    BCPosCoeff______,
    BCVelCoeff______,
    ColorMode_______,
    ColorFactor_____,
    ColorDecay______,
    VisuScale_______,
    VisuSimple______,
    VisuHideOOB_____,
    VerboseLevel____,
  };

  // Particles
  std::vector<Vec::Vec3<float>> Ref;
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<Vec::Vec3<float>> Acc;
  std::vector<Vec::Vec3<float>> For;
  std::vector<Vec::Vec3<float>> Col;
  std::vector<int> Sensor;
  std::vector<int> BCPos;
  std::vector<int> BCVel;
  std::vector<int> BCFor;

  // Buckets for spatial partition
  std::vector<std::vector<std::vector<std::vector<int>>>> Buckets;
  int nX;
  int nY;
  int nZ;

  // Force law
  std::vector<float> ForceLaw;
  float ForceLawStep;

  // Misc
  float SimTime;

  // Scenario functions
  void BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud);
  void BuildScenario(const std::vector<Vec::Vec3<float>>& iPointCloud);

  // Material properties functions
  void BuildForceLaw();

  // Spatial partition functions
  void ComputeBuckets();
  void GetBucketIdx(const Vec::Vec3<float>& iPos, int& oIdxX, int& oIdxY, int& oIdxZ);

  // Simulator functions
  void ComputeForces();
  void ApplyBCForces();
  void StepSimulation();

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
