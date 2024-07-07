#pragma once

// Standard lib
#include <array>
#include <vector>


// Procedural generation of voxel scenes based on Markov algorithm
// - Substitution rule sets generated in code for various scenarios or by importing bitmaps
// - Matching and substitution operated as brute force search on the domain
// - Voxels visualized with shading from hacky volumetric obscuration
//
// Reference on Markov algorithm
// https://en.wikipedia.org/wiki/Markov_algorithm
// https://github.com/mxgmn/MarkovJunior
class MarkovProcGene
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    Scenario________,
    ResolutionX_____,
    ResolutionY_____,
    ResolutionZ_____,
    RuleSizeX_______,
    RuleSizeY_______,
    RuleSizeZ_______,
    NbSubsti________,
    ShadeCoeff______,
    VerboseLevel____,
  };

  int nbX;
  int nbY;
  int nbZ;

  int activeSet;
  int activeRul;

  std::vector<std::vector<std::vector<int>>> Field;
  std::vector<std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>> Dict;

  void FillRuleBox(std::array<std::vector<std::vector<std::vector<int>>>, 2>& ioRule,
                   const int iMinX, const int iMinY, const int iMinZ,
                   const int iMaxX, const int iMaxY, const int iMaxZ,
                   const int iVal, const bool iFillI, const bool iFillO);
  std::array<std::vector<std::vector<std::vector<int>>>, 2> BuildColorSwap(const int iOldColor, const int iNewColor,
                                                                           const std::array<std::vector<std::vector<std::vector<int>>>, 2>& iRule);
  std::array<std::vector<std::vector<std::vector<int>>>, 2> BuildSymmetric(const int iDim1, const int iDim2, const int iDim3,
                                                                           const std::array<std::vector<std::vector<std::vector<int>>>, 2>& iRule);

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  MarkovProcGene();

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
