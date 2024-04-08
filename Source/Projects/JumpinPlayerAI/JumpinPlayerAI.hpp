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
    std::vector<std::vector<int>> Pawns;                              // Flag for presence of pawns (Black= -1, White= +1)
    std::vector<std::vector<std::vector<std::array<int, 2>>>> Moves;  // List of possible moves for each starting position
    int Score;                                                        // Evaluated score of the board
  };

  int nW;                     // Dimensions of the board
  int nH;                     // Dimensions of the board
  BoardState MainBoard;       // Currrent state of the board
  std::array<int, 2> Select;  // Coordinates of selected square

  void UpdateMainBoard();
  void ComputeScore(BoardState &ioBoard);
  void ComputeMoves(BoardState &ioBoard);
  void ComputeDestinations(BoardState &ioBoard,
                           const std::vector<std::vector<bool>> &iOccup,
                           std::vector<std::vector<bool>> &ioVisit,
                           const int iStartW, const int iStartH,
                           const int iW, const int iH);

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
