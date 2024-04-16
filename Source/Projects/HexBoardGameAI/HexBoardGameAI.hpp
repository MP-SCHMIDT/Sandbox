#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Jumpin board game AI
class HexBoardGameAI
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    BoardW__________,
    BoardH__________,
    MoveStreakRed___,
    MoveStreakBlu___,
    RandomPawnInit__,
    ______________00,
    MaxSearchDepth__,
    MaxThinkTime____,
    MaxTreeBoards___,
    MoveSortScore___,
    MoveSortNash____,
    MoveSortRand____,
    ABPruning_______,
    IterDeepening___,
    ______________01,
    ValEdgeConnect__,
    ValCornConnect__,
    ______________02,
    BotStrategyRed__,
    BotStrategyBlu__,
    ______________03,
    ColorFactor_____,
    ______________04,
    TestParamHEX_0__,
    TestParamHEX_1__,
    TestParamHEX_2__,
    TestParamHEX_3__,
    TestParamHEX_4__,
    TestParamHEX_5__,
    ______________05,
    VerboseLevel____,
  };

  // Structure to store a board state
  struct BoardState
  {
    std::vector<BoardState *> SubBoards;  // List of possible boards reachable from the current position sorted in Nash order
    std::vector<std::vector<int>> Pawns;  // Flag grid for presence of pawns on the board (Blu= -1, Red= +1)
    std::array<int, 2> Move;              // Last move that led to this board (or -1 -1 if null move)
    int Score;                            // Evaluated score of the board
    int NashScore;                        // Nash score found in sub tree
    int NashNbSteps;                      // Number of steps to reach the Nash score found in sub tree
  };

  int nW;                                            // Dimensions of the board
  int nH;                                            // Dimensions of the board
  int idxTurn;                                       // Counter for the current turn in the game
  int streakRed;                                     //
  int streakBlu;                                     //
  int turnPeriod;                                    //
  int nbTreeBoards;                                  // Counter for the number of boards evaluated in the current search
  double thinkTime;                                  // Time spent in the last search
  BoardState *RootBoard= nullptr;                    // Current state of the board
  std::vector<std::vector<Vec::Vec3<float>>> Cells;  // Cell centroid coordinates for display and picking
  float cellSize;                                    // Size of a cell
  // Draw
  void DrawBoardTree(const BoardState *iBoard, const int iDepth, const int iDrawMode,
                     const float px, const float py, const float pz,
                     const float radius, const float arcBeg, const float arcEnd);
  void PlotData();

  // Board creation and destruction
  BoardState *CreateBoard(const std::vector<std::vector<int>> &iPawns,
                          const std::array<int, 2> &iMove,
                          const int iDepth);
  void DeleteBoard(BoardState *ioBoard);

  // Evaluation
  bool IsRedTurn(const int iDepth);
  void ComputeBoardScore(BoardState *ioBoard);

  // Search
  void ComputeGameTreeSearch();
  int RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth, int iAlpha, int iBeta);
  void FindPossibleMoves(BoardState *ioBoard, std::vector<std::array<int, 2>> &ioMoves);
  void SortSubBoards(BoardState *ioBoard, const int iDepth, const int iMode);
  int GetIdxBestSubBoard(BoardState *ioBoard, const int iDepth, const int iMode);
  int CompareBoardPair(const BoardState *iBoardA, const BoardState *iBoardB, const int iDepth, const int iMode);

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  HexBoardGameAI();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
