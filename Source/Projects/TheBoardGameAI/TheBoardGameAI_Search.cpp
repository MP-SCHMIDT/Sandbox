#include "TheBoardGameAI.hpp"


// Standard lib

// Algo headers
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void TheBoardGameAI::ComputeGameTreeSearch(const int iMaxDepth) {
  // Start the timer and board counter
  Timer::PushTimer();
  nbTreeBoards= 1;

  // Reset root board
  Field2D<int> Pawns= RootBoard->Pawns;
  DeleteBoard(RootBoard);
  RootBoard= CreateBoard(Pawns, std::vector<std::array<int, 2>>(), 0);
  if (D.UI[GameMode________].I() == 0) ComputeBoardScoreHex(RootBoard);
  if (D.UI[GameMode________].I() == 1) ComputeBoardScoreJmp(RootBoard);
  if (D.UI[GameMode________].I() == 2) ComputeBoardScoreChk(RootBoard);

  // Run the recursive search
  int pruningAlpha= -INT_MAX;
  int pruningBeta= INT_MAX;
  if (D.UI[IterDeepening___].B()) {
    for (int iterMaxSearchDepth= 1; iterMaxSearchDepth <= iMaxDepth; iterMaxSearchDepth++) {
      if ((D.UI[MaxThinkTime____].D() == 0.0 || Timer::CheckTimer() < D.UI[MaxThinkTime____].D()) &&
          (D.UI[MaxTreeBoards___].I() == 0 || nbTreeBoards < D.UI[MaxTreeBoards___].I())) {
        RecursiveTreeSearch(RootBoard, 0, iterMaxSearchDepth, pruningAlpha, pruningBeta);
      }
    }
  }
  else {
    RecursiveTreeSearch(RootBoard, 0, iMaxDepth, pruningAlpha, pruningBeta);
  }

  // Stop the timer
  thinkTime= Timer::PopTimer();
}


int TheBoardGameAI::RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth, int iAlpha, int iBeta) {
  if (ioBoard == nullptr) printf("[ERROR] RecursiveTreeSearch on a null board\n");

  // Handle leaf board
  if (iDepth >= iMaxDepth || ioBoard->Score == INT_MAX || ioBoard->Score == -INT_MAX) {
    ioBoard->NashScore= ioBoard->Score;
    ioBoard->NashNbSteps= 0;
    return ioBoard->NashScore;
  }

  // Create and score the sub boards if none exist
  if (ioBoard->SubBoards.empty()) {
    // Find all the possible moves for the current turn and board
    std::vector<std::vector<std::array<int, 2>>> Moves;
    if (D.UI[GameMode________].I() == 0) FindPossibleMovesHex(ioBoard, iDepth, Moves);
    else if (D.UI[GameMode________].I() == 1) FindPossibleMovesJmp(ioBoard, iDepth, Moves);
    else if (D.UI[GameMode________].I() == 2) FindPossibleMovesChk(ioBoard, iDepth, Moves);

    // Create and score the sub board for each move
    for (std::vector<std::array<int, 2>> move : Moves) {
      BoardState *newBoard= CreateBoard(ioBoard->Pawns, move, iDepth + 1);
      if (D.UI[GameMode________].I() == 0) {
        ExecuteMoveHex(newBoard, iDepth);
        ComputeBoardScoreHex(newBoard);
      }
      if (D.UI[GameMode________].I() == 1) {
        ExecuteMoveJmp(newBoard, iDepth);
        ComputeBoardScoreJmp(newBoard);
      }
      if (D.UI[GameMode________].I() == 2) {
        ExecuteMoveChk(newBoard, iDepth);
        ComputeBoardScoreChk(newBoard);
      }
      ioBoard->SubBoards.push_back(newBoard);
      nbTreeBoards++;
    }

    // Sort the sub boards according to score
    if (D.UI[MoveSortScore___].B()) {
      SortSubBoards(ioBoard, iDepth, 0);
    }
  }
  // Sort the sub boards according to Nash score
  else if (D.UI[MoveSortNash____].B()) {
    SortSubBoards(ioBoard, iDepth, 1);
  }

  // Recursively search the sub boards
  ioBoard->NashScore= IsRedTurn(iDepth) ? -INT_MAX : INT_MAX;
  ioBoard->NashNbSteps= INT_MAX;
  int prunedMaxDepth= iMaxDepth;
  for (int k= 0; k < (int)ioBoard->SubBoards.size(); k++) {
    const int score= RecursiveTreeSearch(ioBoard->SubBoards[k], iDepth + 1, prunedMaxDepth, iAlpha, iBeta);
    if ((IsRedTurn(iDepth) && ioBoard->NashScore < score) ||
        (!IsRedTurn(iDepth) && ioBoard->NashScore > score) ||
        (ioBoard->NashScore == score && ioBoard->NashNbSteps > ioBoard->SubBoards[k]->NashNbSteps + 1)) {
      ioBoard->NashScore= score;
      ioBoard->NashNbSteps= ioBoard->SubBoards[k]->NashNbSteps + 1;
    }
    if (D.UI[ABPruning_______].B()) {
      if (IsRedTurn(iDepth)) {
        iAlpha= std::max(iAlpha, ioBoard->NashScore);
        if (ioBoard->NashScore >= iBeta) prunedMaxDepth= 0;
      }
      else {
        iBeta= std::min(iBeta, ioBoard->NashScore);
        if (ioBoard->NashScore <= iAlpha) prunedMaxDepth= 0;
      }
    }
  }

  return ioBoard->NashScore;
}


// Sort the sub boards for the picked comparison mode
void TheBoardGameAI::SortSubBoards(BoardState *ioBoard, const int iDepth, const int iMode) {
  if (ioBoard == nullptr) printf("[ERROR] SortSubBoards on a null board\n");

  // Bubble sort of the sub boards
  int kEnd= (int)ioBoard->SubBoards.size();
  do {
    int kNew= 0;
    for (int k= 1; k < kEnd; k++) {
      const int comparison= CompareBoardPair(ioBoard->SubBoards[k - 1], ioBoard->SubBoards[k], iDepth, iMode);
      if ((comparison < 0) || (comparison == 0 && D.UI[MoveSortRand____].B() && Random::Val(0, 1) == 0)) {
        BoardState *tmp= ioBoard->SubBoards[k - 1];
        ioBoard->SubBoards[k - 1]= ioBoard->SubBoards[k];
        ioBoard->SubBoards[k]= tmp;
        kNew= k;
      }
    }
    kEnd= kNew;
  } while (kEnd > 1);
}


// Get the index of the best sub board for the picked comparison mode
int TheBoardGameAI::GetIdxBestSubBoard(BoardState *ioBoard, const int iDepth, const int iMode) {
  if (ioBoard == nullptr) printf("[ERROR] GetIdxBestSubBoard on a null board\n");

  // Get the index of the best sub board
  int idxBest= 0;
  for (int k= 1; k < (int)ioBoard->SubBoards.size(); k++) {
    const int comparison= CompareBoardPair(ioBoard->SubBoards[idxBest], ioBoard->SubBoards[k], iDepth, iMode);
    if ((comparison < 0) || (comparison == 0 && D.UI[MoveSortRand____].B() && Random::Val(0, 1) == 0))
      idxBest= k;
  }
  return idxBest;
}


// Compare a board pair for the picked comparison mode
// A best= +1, B best= -1, Equal= 0
int TheBoardGameAI::CompareBoardPair(const BoardState *iBoardA, const BoardState *iBoardB, const int iDepth, const int iMode) {
  if (iBoardA == nullptr || iBoardB == nullptr) printf("[ERROR] CompareBoardPair on a null board\n");

  if (iMode == 0) {
    if ((IsRedTurn(iDepth) && iBoardA->Score > iBoardB->Score) ||
        (!IsRedTurn(iDepth) && iBoardA->Score < iBoardB->Score))
      return +1;
    if (iBoardA->Score == iBoardB->Score)
      return 0;
  }
  else {
    if ((IsRedTurn(iDepth) && iBoardA->NashScore > iBoardB->NashScore) ||
        (!IsRedTurn(iDepth) && iBoardA->NashScore < iBoardB->NashScore) ||
        (iBoardA->NashScore == iBoardB->NashScore && iBoardA->NashNbSteps < iBoardB->NashNbSteps))
      return +1;
    if (iBoardA->NashScore == iBoardB->NashScore && iBoardA->NashNbSteps == iBoardB->NashNbSteps)
      return 0;
  }
  return -1;
}
