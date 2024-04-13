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
    BotBluPlayer____,
    BotRedPlayer____,
    StartingRows____,
    ______________00,
    MaxSearchDepth__,
    MaxThinkTime____,
    MaxTreeNodes____,
    ______________01,
    ValPushTotal____,
    ValPushLast_____,
    ValSoftStranded_,
    ValHardStranded_,
    ______________02,
    ColorFactor_____,
    ______________03,
    TestParamGAI_0__,
    TestParamGAI_1__,
    TestParamGAI_2__,
    TestParamGAI_3__,
    TestParamGAI_4__,
    TestParamGAI_5__,
    ______________04,
    VerboseLevel____,
  };

  // Structure to store a board state
  struct BoardState
  {
    std::vector<BoardState *> SubBoards;  // List of possible boards reachable from the current position sorted in Nash order
    std::vector<std::vector<int>> Pawns;  // Flag grid for presence of pawns on the board (Red= -1, Blu= +1)
    std::array<int, 4> Move;              // Move source and destination that leads to this board (or 0 0 0 0 if root board)
    int Score;                            // Evaluated score of the board
    int NashScore;                        // Nash score found in sub tree
    int NashNbSteps;                      // Number of steps to reach the Nash score found in sub tree
  };

  int nW;                          // Dimensions of the board
  int nH;                          // Dimensions of the board
  int idxTurn;                     // Counter for the current turn in the game
  int nbTreeNodes;                 // Counter for the number of nodes evaluated in the current search
  double thinkTime;                // Time spent in the last search
  BoardState *RootBoard= nullptr;  // Current state of the board
  int wSel;                        // Coordinates of the currently selected pawn
  int hSel;                        // Coordinates of the currently selected pawn

  // Draw
  void DrawBoardTree(const BoardState *iBoard, const int iDepth,
                     const float px, const float py, const float pz,
                     const float radius, const float arcBeg, const float arcEnd);
  void PlotData();

  // Board creation and destruction
  BoardState *CreateBoard(const std::vector<std::vector<int>> &iPawns,
                          const std::array<int, 4> &iMove);
  void DeleteBoard(BoardState *ioBoard);

  // Evaluation
  void ComputeBoardScore(BoardState *ioBoard);

  // Search
  void ComputeGameTreeSearch();
  void RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth);
  void RecursivePawnJumps(BoardState *ioBoard, const int iDepth,
                          const int iStartW, const int iStartH,
                          const int iCurrW, const int iCurrH,
                          std::vector<std::vector<bool>> &ioVisit,
                          std::vector<std::array<int, 2>> &ioJumps);

  // Utility
  void inline SortedBoardInsertion(BoardState *ioBoard, const int iDepth, BoardState *NewBoard);
  bool inline IsBluTurn(const int iDepth);

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
