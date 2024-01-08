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
    ForceColli______,
    ForcePartic_____,
    ForceAdhesion___,
    ForceBoundary___,
    ForceDamping____,
    ForceGravity____,
    ForceAdhMode____,
    ForceAdhCoeff___,
    ColorMode_______,
    VerboseLevel____,
  };

  // Rules
  std::vector<double> ParticleTypeMass;
  std::vector<double> ParticleTypeRadius;
  std::vector<std::vector<double>> ParticleTypeAmpli;
  std::vector<std::vector<double>> ParticleTypeReach;

  // Particles
  std::vector<int> ParticleType;
  std::vector<Vec::Vec3<double>> ParticlePos;
  std::vector<Vec::Vec3<double>> ParticleVel;
  std::vector<Vec::Vec3<double>> ParticleAcc;
  std::vector<Vec::Vec3<double>> ParticleFor;

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
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
