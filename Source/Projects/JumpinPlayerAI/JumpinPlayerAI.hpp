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
    StartingRows____,
    VerboseLevel____,
  };

  // Structure to store a board state
  struct BoardState
  {
    std::vector<std::vector<char>> Pawns;                              // Flag for presence of pawns (Black= -1, White= +1)
    std::vector<std::vector<std::vector<std::array<char, 2>>>> Moves;  // List of possible moves for each starting position
    int Score;                                                         // Evaluated score of the board
  };

  int nW;                      // Dimensions of the board
  int nH;                      // Dimensions of the board
  BoardState Board;            // Currrent state of the board
  std::array<char, 2> Select;  // Coordinates of selected square

  void ComputeScore();
  void ComputeMoves();
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
