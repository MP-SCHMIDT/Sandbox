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
    EnableBluPlayer_,
    EnableRedPlayer_,
    StartingRows____,
    ______________00,
    SearchDepth_____,
    ______________01,
    ColorFactor_____,
    ______________02,
    VerboseLevel____,
  };

  // Structure to store a board state
  struct BoardState
  {
    std::vector<std::vector<int>> Pawns;                              // Flag for presence of pawns (Red= -1, Blu= +1)
    std::vector<std::vector<std::vector<std::array<int, 2>>>> Moves;  // List of possible moves for each starting position
    std::vector<std::vector<std::vector<BoardState *>>> SubBoards;    // Resulting board states for each potential move
    int Score;                                                        // Evaluated score of the board
    int ScoreBest;                                                    // Best score found in sub tree
    int StepsBest;                                                    // Number of steps to reach the best score found in sub tree
    int WinState;                                                     // Number indicating if a player has a guaranteed win
    bool BestIsSet;                                                   // Flag indicating if the best move is set
  };

  int nW;                          // Dimensions of the board
  int nH;                          // Dimensions of the board
  int idxTurn;                     // Counter for the current turn in the game
  BoardState *RootBoard= nullptr;  // Currrent state of the board
  int wSel;                        // Coordinates of the currently selected pawn
  int hSel;                        // Coordinates of the currently selected pawn

  // Draw
  void DrawBoardTree(const BoardState *iBoard, const int iDepth,
                     const float px, const float py, const float pz,
                     const float radius, const float arcBeg, const float arcEnd);
  void PlotData();

  // Board creation and destruction
  BoardState *CreateBoard();
  void DeleteBoard(BoardState *ioBoard);

  // Game tree board computation
  void ComputeGameTreeSearch(BoardState *ioBoard, const int iDepth);
  void CheckWinState(BoardState *ioBoard, const int iDepth);
  void ComputeBoardScore(BoardState *ioBoard);
  void ComputePawnJumps(BoardState *ioBoard,
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
