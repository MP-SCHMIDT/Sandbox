#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Vec.hpp"


// Simulation with a Position-Based Dynamics (PBD) approach
// - Node positions are corrected by solving position constraints
// - Velocities are deduced from changed positions
// - Particle positions are iteratively updated with explicit time integration
// - Collision constraints are resolved by correcting positions of particles with a Gauss Seidel relaxation
// - Masses (and inertia) are lumped into nodes
// - Constraints: pairwise node collisions, edge length, triangle area, tetrahedron volume, mouse repulsion
//
// References for subset of functionalities
// - Tutorials 9 & 10 from Matthias Müller
// - https://matthias-research.github.io/pages/tenMinutePhysics/index.html
// - https://www.youtube.com/watch?v=jrociOAYqxA
// - https://www.youtube.com/watch?v=lS_qeBy3aQI
class PosiBasedDynam
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    DomainX_________,
    DomainY_________,
    DomainZ_________,
    NumParticl______,
    RadParticl______,
    RadMouse________,
    AddTetMeshID____,
    ______________00,
    TimeStep________,
    NbSubSteps______,
    MaterialDensity_,
    ForceDrag_______,
    ForceGrav_______,
    StiffBox________,
    StiffCollision__,
    StiffEdgeLength_,
    StiffTriArea____,
    StiffTetVolume__,
    StiffMouseBall__,
    RandPerturb_____,
    ______________01,
    ColorMode_______,
    ColorFactor_____,
    VisuSimple______,
    ______________02,
    TestParamPBD_0__,
    TestParamPBD_1__,
    TestParamPBD_2__,
    TestParamPBD_3__,
    TestParamPBD_4__,
    VerboseLevel____,
  };

  std::vector<Vec::Vec3<float>> PosOld;
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<std::array<int, 2>> Edge;
  std::vector<std::array<int, 3>> Tri;
  std::vector<std::array<int, 4>> Tet;
  std::vector<float> EdgeLen;
  std::vector<float> TriArea;
  std::vector<float> TetVol;
  std::vector<float> Mass;
  std::vector<float> MassInv;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  PosiBasedDynam();

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

  void TimeIntegrate();

  void ApplyBoxDomainConstraint(const float iStiffness, const float iTimestep);
  void ApplyNodeCollisionConstraint(const float iStiffness, const float iTimestep);
  void ApplyEdgeLengthConstraint(const float iStiffness, const float iTimestep);
  void ApplyTriangleAreaConstraint(const float iStiffness, const float iTimestep);
  void ApplyTetrahedronVolumeConstraint(const float iStiffness, const float iTimestep);
  void ApplyMouseBallConstraint(const float iStiffness, const float iTimestep);
  
  void DrawScene();

  // Utility functions
  inline float GetTetVolume(const Vec::Vec3<float> &x0, const Vec::Vec3<float> &x1, const Vec::Vec3<float> &x2, const Vec::Vec3<float> &x3) {
    return (((x1-x0).cross(x2-x0)).dot(x3-x0)) / 6.0f;
  }
  inline float GetTriArea(const Vec::Vec3<float> &x0, const Vec::Vec3<float> &x1, const Vec::Vec3<float> &x2) {
    return 0.5f * ((x1-x0).cross(x2-x0)).norm();
  }
  inline float GetEdgeLength(const Vec::Vec3<float> &x0, const Vec::Vec3<float> &x1) {
    return (x1-x0).norm();
  }
};
