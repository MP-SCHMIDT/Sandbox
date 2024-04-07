#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Jumpin board game AI
class JumpinPlayerAI
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    BoardW__________,
    BoardH__________,
    VerboseLevel____,
  };

  int nW;                                       // Dimensions of the board
  int nH;                                       // Dimensions of the board
  std::vector<std::vector<bool>> White;         // Flag for presence of white pieces
  std::vector<std::vector<bool>> Black;         // Flag for presence of black pieces
  std::vector<std::vector<bool>> Occupied;      // Flag for presence of any pieces
  std::vector<std::vector<bool>> Destinations;  // Flag for possible destinations of picked square
  std::array<int, 2> Select;                    // Coordinates of selected square

  void ComputeDestinations(const int w, const int h);

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  JumpinPlayerAI();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
