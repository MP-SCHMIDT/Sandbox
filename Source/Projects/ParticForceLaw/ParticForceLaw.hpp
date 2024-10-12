#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Type/Field.hpp"
#include "Type/Vec.hpp"


// Particle-based physics simulation
// - Memoryless Isotropic Point Particles defined only by position and velocity
// - Different from most DEM because no internal variables, history, rotations, etc.
// - Base particle cloud as FCC lattice, Poisson sphere sampling, etc.
// - Isotropic force law defining force vs distance (hard-coded presets or UI params)
// - Explicit time integration with forward Euler
// - Boundary conditions as simple overwrite or as additional force term
// - O(n) neighborhood search with spatial partition into buckets of fixed size
// - Multimaterial interaction by average of forces laws at material interfaces
// - Isosurface visualization with metaball scalar field and marching cubes
// - Expected performance: ~1M-2M particle updates per second on laptop CPU (Buckets 32k x 40)
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
    ConstrainDim2D__,
    ______________00,
    StepsPerDraw____,
    TimeStep________,
    BucketCapacity__,
    BucketFillCoeff_,
    DampingRadRel___,
    DampingVelRel___,
    MaterialDensity_,
    ______________01,
    BCVelX__________,
    BCVelY__________,
    BCVelZ__________,
    BCForX__________,
    BCForY__________,
    BCForZ__________,
    UseForceControl_,
    BCPosCoeff______,
    BCVelCoeff______,
    ______________02,
    ForceLawPresetA_,
    ForceLawPresetB_,
    ForceLawScaleA__,
    ForceLawScaleB__,
    ______________03,
    ForceLawA_0_00__,
    ForceLawA_0_80__,
    ForceLawA_0_90__,
    ForceLawA_0_95__,
    ForceLawA_1_00__,
    ForceLawA_1_05__,
    ForceLawA_1_10__,
    ForceLawA_1_20__,
    ForceLawA_1_30__,
    ForceLawA_1_40__,
    ForceLawA_1_50__,
    ForceLawA_2_00__,
    ForceLawA_2_50__,
    ForceLawA_3_00__,
    ______________04,
    ForceLawB_0_00__,
    ForceLawB_0_80__,
    ForceLawB_0_90__,
    ForceLawB_0_95__,
    ForceLawB_1_00__,
    ForceLawB_1_05__,
    ForceLawB_1_10__,
    ForceLawB_1_20__,
    ForceLawB_1_30__,
    ForceLawB_1_40__,
    ForceLawB_1_50__,
    ForceLawB_2_00__,
    ForceLawB_2_50__,
    ForceLawB_3_00__,
    ______________05,
    MetaballVoxSize_,
    MetaballIsoval__,
    ColorMode_______,
    ColorFactor_____,
    VisuScale_______,
    VisuSimple______,
    VisuHideOOB_____,
    VisuMinNeighbor_,
    ______________06,
    TestParamMIP_0__,
    TestParamMIP_1__,
    TestParamMIP_2__,
    TestParamMIP_3__,
    TestParamMIP_4__,
    TestParamMIP_5__,
    VerboseLevel____,
  };

  // Particles
  std::vector<Vec::Vec3<float>> Ref;
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<Vec::Vec3<float>> For;
  std::vector<Vec::Vec3<float>> Col;
  std::vector<float> ForceMag;
  std::vector<int> Neighbors;
  std::vector<int> Mat;
  std::vector<int> Sensor;
  std::vector<int> BCPos;
  std::vector<int> BCVel;
  std::vector<int> BCFor;

  // Buckets for spatial partition
  Field::Field3<std::vector<int>> Buckets;
  bool BucketOverflown;

  // Metaball visualization
  bool MetaballIsUpdated;
  std::vector<std::array<double, 3>> Verts;
  std::vector<std::array<int, 3>> Tris;

  // Force law
  std::vector<std::vector<float>> ForceLaw;
  std::vector<float> ForceLawStep;
  std::vector<float> ForceLawRange;

  // Misc
  float SimTime;
  int RunID;

  // Scenario functions
  void BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud);
  void BuildScenario(const std::vector<Vec::Vec3<float>>& iPointCloud);

  // Material properties functions
  void BuildForceLaws();
  void BuildForceLawPolyline(const double v0_00, const double v0_80, const double v0_90, const double v0_95,
                             const double v1_00, const double v1_05, const double v1_10, const double v1_20,
                             const double v1_30, const double v1_40, const double v1_50, const double v2_00,
                             const double v2_50, const double v3_00, const int iIdxMat);

  // Spatial partition functions
  void ComputeBuckets();
  void GetBucketIdx(const Vec::Vec3<float>& iPos, int& oIdxX, int& oIdxY, int& oIdxZ);

  // Simulator functions
  void ComputeForces();
  void ApplyBCForces();
  void StepSimulation();

  // Visualization functions
  void ComputeMetaballs();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  ParticForceLaw();

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
