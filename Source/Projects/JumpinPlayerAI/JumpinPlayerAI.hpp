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
    SinglePlayer____,
    SearchDepth_____,
    TreeStepDist____,
    TreeStepRadians_,
    TreeFactDist____,
    TreeFactRadians_,
    ColorFactor_____,
    VerboseLevel____,
  };

  // Structure to store a board state
  struct BoardState
  {
    std::vector<std::vector<int>> Pawns;                              // Flag for presence of pawns (Black= -1, White= +1)
    std::vector<std::vector<std::vector<std::array<int, 2>>>> Moves;  // List of possible moves for each starting position
    std::vector<std::vector<std::vector<BoardState *>>> SubBoards;    // Resulting board states for each potential move
    int Score;                                                        // Evaluated score of the board
    int ScoreBestWhite;                                               // Best score for white found in sub tree
    int ScoreBestBlack;                                               // Best score for black found in sub tree
    int StepsBestWhite;                                               // Number of steps to reach best score for white found in sub tree
    int StepsBestBlack;                                               // Number of steps to reach best score for black found in sub tree
    std::array<std::array<int, 2>, 2> MoveBestWhite;                  // Move for the best score for white found in sub tree
    std::array<std::array<int, 2>, 2> MoveBestBlack;                  // Move for the best score for black found in sub tree
  };

  int nW;                          // Dimensions of the board
  int nH;                          // Dimensions of the board
  int idxTurn;                     // Counter for the current turn in the game
  BoardState *RootBoard= nullptr;  // Currrent state of the board
  int wSel;                        // Coordinates of the currently selected pawn
  int hSel;                        // Coordinates of the currently selected pawn

  // Draw
  void DrawBoardTree(const BoardState *iBoard, const int iDrawMode,
                     const float px, const float py, const float pz,
                     const float dist, const float radians,
                     const float distStep, const float radiansStep,
                     const float distFact, const float radiansFact);

  // Board creation and destruction
  BoardState *CreateBoard();
  void DeleteBoard(BoardState *ioBoard);
  void DeleteSubBoards(BoardState *ioBoard, const int w, const int h);

  // Board computation
  void ComputeBoardScore(BoardState *ioBoard);
  void ConsolidateBoardScores(BoardState *ioBoard, const int iDepth);
  void ComputeBoardMoves(BoardState *ioBoard, const int iDepth);
  void ComputePawnDestinations(BoardState *ioBoard,
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
