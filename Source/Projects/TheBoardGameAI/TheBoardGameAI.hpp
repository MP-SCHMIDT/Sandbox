#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Jumpin board game AI
class TheBoardGameAI
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    GameMode________,
    BoardW__________,
    BoardH__________,
    MoveStreakRed___,
    MoveStreakBlu___,
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
    HexEdgeConnect__,
    HexCornConnect__,
    JmpPushTotal____,
    JmpPushLast_____,
    JmpSoftStranded_,
    JmpHardStranded_,
    ChkMaterial_____,
    ______________02,
    RandomMoves_____,
    BotStrategyRed__,
    BotStrategyBlu__,
    ______________03,
    ColorMode_______,
    ColorFactor_____,
    ______________04,
    TestParamGAI_0__,
    TestParamGAI_1__,
    TestParamGAI_2__,
    TestParamGAI_3__,
    TestParamGAI_4__,
    TestParamGAI_5__,
    ______________05,
    VerboseLevel____,
  };

  // Structure to store a board state
  struct BoardState
  {
    std::vector<BoardState *> SubBoards;   // List of possible boards reachable from the current position sorted in Nash order
    std::vector<std::vector<int>> Pawns;   // Flag grid for presence of pawns on the board (Blu= -1, Red= +1)
    std::vector<std::array<int, 2>> Move;  // Last move sequence that led to this board
    int Score;                             // Evaluated score of the board
    int NashScore;                         // Nash score found in sub tree
    int NashNbSteps;                       // Number of steps to reach the Nash score found in sub tree
  };

  int nW;                          // Dimensions of the board
  int nH;                          // Dimensions of the board
  int idxTurn;                     // Counter for the current turn in the game
  int streakRed;                   // Number of consecutive Red moves chosen for the game sequencing
  int streakBlu;                   // Number of consecutive Blu moves chosen for the game sequencing
  int nbTreeBoards;                // Counter for the number of boards evaluated in the current search
  double thinkTime;                // Time spent in the last search
  BoardState *RootBoard= nullptr;  // Current state of the board
  int wSel;                        // Coordinates of the currently selected pawn
  int hSel;                        // Coordinates of the currently selected pawn

  // Hex draw data
  std::vector<std::vector<Vec::Vec3<float>>> Cells;  // Cell centroid coordinates for display and picking
  float cellSize;                                    // Size of a cell

  // Draw
  void DrawBoardHex();
  void DrawBoardJmp();
  void DrawBoardChk();
  void DrawBoardTree(const BoardState *iBoard, const int iDepth, const int iDrawMode,
                     const float px, const float py, const float pz,
                     const float radius, const float arcBeg, const float arcEnd);
  void PlotData();

  // Board creation and destruction
  BoardState *CreateBoard(const std::vector<std::vector<int>> &iPawns, const int iDepth);
  void DeleteBoard(BoardState *ioBoard);

  // Evaluation
  bool IsRedTurn(const int iDepth);
  void ComputeBoardScoreHex(BoardState *ioBoard);
  void ComputeBoardScoreJmp(BoardState *ioBoard);
  void ComputeBoardScoreChk(BoardState *ioBoard);

  // Search
  void ComputeGameTreeSearch(const int iMaxDepth);
  int RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth, int iAlpha, int iBeta);
  void SortSubBoards(BoardState *ioBoard, const int iDepth, const int iMode);
  int GetIdxBestSubBoard(BoardState *ioBoard, const int iDepth, const int iMode);
  int CompareBoardPair(const BoardState *iBoardA, const BoardState *iBoardB, const int iDepth, const int iMode);

  // Move generation
  void FindPossibleMovesHex(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves);
  void FindPossibleMovesJmp(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves);
  void FindPossibleMovesChk(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves);
  void JmpRecursivePawnMoves(BoardState *ioBoard,
                             const int iJumpW, const int iJumpH,
                             std::vector<std::vector<bool>> &ioVisit,
                             std::vector<std::array<int, 2>> &ioDestinations);

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  TheBoardGameAI();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void Refresh();
  void KeyPress(const unsigned char key);
  void MousePress(const unsigned char mouse);
  void Animate();
  void Draw();
};