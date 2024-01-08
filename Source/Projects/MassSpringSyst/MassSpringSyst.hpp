#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Mass-spring soft body simulation
// - TODO
//
// Reference mass spring
// https://www.cs.rpi.edu/~cutler/classes/advancedgraphics/S17/lectures/06_mass_spring_systems.pdf
//
// Reference implicit solve
// https://njoubert.com/assets/mine_simulation.pdf
// http://www.cs.cmu.edu/~baraff/sigcourse/notese.pdf
// http://www.cs.cmu.edu/~baraff/sigcourse/slidese.pdf
// https://blog.mmacklin.com/2012/05/04/implicitsprings/
// https://stackoverflow.com/questions/3897424/implementing-semi-implicit-backward-euler-in-a-1-dof-mass-spring-system?rq=4
// https://hugi.scene.org/online/hugi28/hugi%2028%20-%20coding%20corner%20uttumuttu%20implementing%20the%20implicit%20euler%20method%20for%20mass-spring%20systems.htm
class MassSpringSyst
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    Scenario________,
    InputFile_______,
    DomainX_________,
    DomainY_________,
    DomainZ_________,
    DistribMode_____,
    NbNodesTarg_____,
    LinkDist________,
    TimeStep________,
    IntegMode_______,
    SolvMaxIter_____,
    CoeffExt________,
    CoeffGravi______,
    CoeffSpring_____,
    CoeffDamp_______,
    ColorFactor_____,
    ColorMode_______,
    VerboseLevel____,
  };

  int N;
  std::vector<std::vector<int>> Adj;
  std::vector<Vec::Vec3<float>> Ref;
  std::vector<Vec::Vec3<float>> Pos;
  std::vector<Vec::Vec3<float>> Vel;
  std::vector<Vec::Vec3<float>> Acc;
  std::vector<Vec::Vec3<float>> For;
  std::vector<Vec::Vec3<float>> Ext;
  std::vector<Vec::Vec3<float>> Fix;
  std::vector<float> Mas;

  void ComputeForces();
  void StepForwardInTime();
  void ApplyBCPos();
  void ApplyBCVel();
  void ApplyBCFor();

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  MassSpringSyst();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
