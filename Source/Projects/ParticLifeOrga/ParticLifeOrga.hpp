#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Particle Life organization
// TODO explain
//
// TODO link PL video
class ParticLifeOrga
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    NbTypes_________,
    NbPartic________,
    TimeStep________,
    NumIntegMode____,
    Use2DProj_______,
    DomainX_________,
    DomainY_________,
    DomainZ_________,
    RuleMinMass_____,
    RuleMaxMass_____,
    RuleMinRadius___,
    RuleMaxRadius___,
    RuleMinAmpli____,
    RuleMaxAmpli____,
    RuleMinReach____,
    RuleMaxReach____,
    RuleAdhReach____,
    ModeAdhesion____,
    ForceColli______,
    ForcePartic_____,
    ForceBoundary___,
    ForceDamping____,
    ForceGravity____,
    ForceAdhesion___,
    ColorMode_______,
    VerboseLevel____,
  };

  // Rules
  std::vector<double> TypeMass;
  std::vector<double> TypeRadius;
  std::vector<std::vector<double>> TypeAmpli;
  std::vector<std::vector<double>> TypeReach;

  // Particles
  std::vector<int> Type;
  std::vector<Vec::Vec3<double>> Pos;
  std::vector<Vec::Vec3<double>> Vel;
  std::vector<Vec::Vec3<double>> Acc;
  std::vector<Vec::Vec3<double>> For;

  void GenerateParticleRules();
  void GenerateParticleCloud();
  void ComputeForces();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  ParticLifeOrga();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void Refresh();
  void KeyPress();
  void MousePress();
  void Animate();
  void Draw();
};
