#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers


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
    RandomBoard_____,
    SinglePlayer____,
    BotStrategyRed__,
    BotStrategyBlu__,
    ______________00,
    MaxSearchDepth__,
    MaxThinkTime____,
    MaxTreeBoards___,
    MoveOrdering____,
    ABPruning_______,
    IterDeepening___,
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
    std::vector<std::vector<int>> Pawns;  // Flag grid for presence of pawns on the board (Blu= -1, Red= +1)
    std::array<int, 4> Move;              // Move source and destination that leads to this board (or 0 0 0 0 if root board)
    int Score;                            // Evaluated score of the board
    int NashScore;                        // Nash score found in sub tree
    int NashNbSteps;                      // Number of steps to reach the Nash score found in sub tree
  };

  int nW;                          // Dimensions of the board
  int nH;                          // Dimensions of the board
  int idxTurn;                     // Counter for the current turn in the game
  int nbTreeBoards;                // Counter for the number of boards evaluated in the current search
  double thinkTime;                // Time spent in the last search
  BoardState *RootBoard= nullptr;  // Current state of the board
  int wSel;                        // Coordinates of the currently selected pawn
  int hSel;                        // Coordinates of the currently selected pawn

  // Draw
  void DrawBoardTree(const BoardState *iBoard, const int iDepth, const int iDrawMode,
                     const float px, const float py, const float pz,
                     const float radius, const float arcBeg, const float arcEnd);
  void PlotData();

  // Board creation and destruction
  BoardState *CreateBoard(const std::vector<std::vector<int>> &iPawns,
                          const std::array<int, 4> &iMove,
                          const int iDepth);
  void DeleteBoard(BoardState *ioBoard);

  // Evaluation
  bool IsRedTurn(const int iDepth);
  void ComputeBoardScore(BoardState *ioBoard);

  // Search
  void ComputeGameTreeSearch();
  int RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth, int &alpha, int &beta);
  void RecursivePawnMoves(BoardState *ioBoard,
                          const int iPawnW, const int iPawnH,
                          const int iJumpW, const int iJumpH,
                          std::vector<std::vector<bool>> &ioVisit,
                          std::vector<std::array<int, 4>> &ioMoves);
  int GetIdxNashSubBoard(BoardState *ioBoard, const int iDepth);
  void SortSubBoards(BoardState *ioBoard, const int iDepth);

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
