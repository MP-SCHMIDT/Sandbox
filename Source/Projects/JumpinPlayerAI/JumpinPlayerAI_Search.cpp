#include "JumpinPlayerAI.hpp"


// Standard lib

// Algo headers
#include "Math/Field.hpp"
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void JumpinPlayerAI::ComputeGameTreeSearch() {
  // Start the timer and board counter
  Timer::PushTimer();
  nbTreeBoards= 1;

  // Reset root board
  std::vector<std::vector<int>> Pawns= RootBoard->Pawns;
  DeleteBoard(RootBoard);
  RootBoard= CreateBoard(Pawns, {-1, -1, -1, -1}, 0);
  ComputeBoardScore(RootBoard);

  // Run the recursive search
  int pruningAlpha= -INT_MAX;
  int pruningBeta= INT_MAX;
  if (D.UI[IterDeepening___].B()) {
    for (int iterMaxSearchDepth= 1; iterMaxSearchDepth <= D.UI[MaxSearchDepth__].I(); iterMaxSearchDepth++) {
      RecursiveTreeSearch(RootBoard, 0, iterMaxSearchDepth, pruningAlpha, pruningBeta);
    }
  }
  else {
    RecursiveTreeSearch(RootBoard, 0, D.UI[MaxSearchDepth__].I(), pruningAlpha, pruningBeta);
  }

  // Stop the timer
  thinkTime= Timer::PopTimer();
}


int JumpinPlayerAI::RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth, int iAlpha, int iBeta) {
  if (ioBoard == nullptr) printf("[ERROR] RecursiveTreeSearch on a null board\n");

  // Handle stopping criterion
  if ((D.UI[MaxThinkTime____].D() > 0.0 && Timer::CheckTimer() >= D.UI[MaxThinkTime____].D()) ||
      (D.UI[MaxTreeBoards___].I() > 0 && nbTreeBoards >= D.UI[MaxTreeBoards___].I())) {
    ioBoard->NashScore= IsRedTurn(iDepth) ? -INT_MAX : INT_MAX;
    ioBoard->NashNbSteps= INT_MAX;
    return ioBoard->NashScore;
  }

  // Handle leaf board
  if (iDepth >= iMaxDepth || ioBoard->Score == INT_MAX || ioBoard->Score == -INT_MAX) {
    ioBoard->NashScore= ioBoard->Score;
    ioBoard->NashNbSteps= 0;
    return ioBoard->NashScore;
  }

  // Create and score the sub boards if none exist
  if (ioBoard->SubBoards.empty()) {
    // Find all the possible pawn moves from the current turn
    std::vector<std::array<int, 4>> Moves;
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if ((ioBoard->Pawns[w][h] > 0 && IsRedTurn(iDepth)) || (ioBoard->Pawns[w][h] < 0 && !IsRedTurn(iDepth))) {
          // Compute possible moves for the current pawn
          std::vector<std::vector<bool>> Visit= Field::AllocField2D(nW, nH, false);
          Visit[w][h]= true;
          const int oldPawn= ioBoard->Pawns[w][h];
          ioBoard->Pawns[w][h]= 0;
          RecursivePawnMoves(ioBoard, w, h, w, h, Visit, Moves);
          ioBoard->Pawns[w][h]= oldPawn;
        }
      }
    }

    // Add the turn skip move if no other possible move exists
    if (Moves.empty())
      Moves.push_back(std::array<int, 4>{-1, -1, -1, -1});

    // Create and score the sub board for each move
    for (std::array<int, 4> move : Moves) {
      BoardState *newBoard= CreateBoard(ioBoard->Pawns, move, iDepth + 1);
      if (move[0] != -1 && move[1] != -1 && move[2] != -1 && move[3] != -1) {
        newBoard->Pawns[move[2]][move[3]]= newBoard->Pawns[move[0]][move[1]];
        newBoard->Pawns[move[0]][move[1]]= 0;
      }
      ComputeBoardScore(newBoard);
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
  for (int k= 0; k < (int)ioBoard->SubBoards.size(); k++) {
    const int score= RecursiveTreeSearch(ioBoard->SubBoards[k], iDepth + 1, iMaxDepth, iAlpha, iBeta);
    if ((IsRedTurn(iDepth) && ioBoard->NashScore < score) ||
        (!IsRedTurn(iDepth) && ioBoard->NashScore > score) ||
        (ioBoard->NashScore == score && ioBoard->NashNbSteps > ioBoard->SubBoards[k]->NashNbSteps + 1)) {
      ioBoard->NashScore= score;
      ioBoard->NashNbSteps= ioBoard->SubBoards[k]->NashNbSteps + 1;
    }
    if (D.UI[ABPruning_______].B()) {
      if (IsRedTurn(iDepth)) {
        iAlpha= std::max(iAlpha, ioBoard->NashScore);
        if (ioBoard->NashScore >= iBeta) break;
      }
      else {
        iBeta= std::min(iBeta, ioBoard->NashScore);
        if (ioBoard->NashScore <= iAlpha) break;
      }
    }
  }

  return ioBoard->NashScore;
}


void JumpinPlayerAI::RecursivePawnMoves(BoardState *ioBoard,
                                        const int iPawnW, const int iPawnH,
                                        const int iJumpW, const int iJumpH,
                                        std::vector<std::vector<bool>> &ioVisit,
                                        std::vector<std::array<int, 4>> &ioMoves) {
  if (ioBoard == nullptr) printf("[ERROR] RecursivePawnMoves on a null board\n");

  for (int idxDir= 0; idxDir < 4; idxDir++) {
    const int wInc= (idxDir == 0) ? (-1) : ((idxDir == 1) ? (+1) : (0));
    const int hInc= (idxDir == 2) ? (-1) : ((idxDir == 3) ? (+1) : (0));
    if (iJumpW + 2 * wInc < 0 || iJumpW + 2 * wInc >= nW ||
        iJumpH + 2 * hInc < 0 || iJumpH + 2 * hInc >= nH) continue;
    if (ioBoard->Pawns[iJumpW + wInc][iJumpH + hInc] == 0) continue;
    int wOff= iJumpW + 2 * wInc;
    int hOff= iJumpH + 2 * hInc;
    while (wOff >= 0 && wOff < nW && hOff >= 0 && hOff < nH) {
      if (ioBoard->Pawns[wOff][hOff] == 0) {
        if (!ioVisit[wOff][hOff]) {
          ioVisit[wOff][hOff]= true;
          ioMoves.push_back(std::array<int, 4>{iPawnW, iPawnH, wOff, hOff});
          RecursivePawnMoves(ioBoard, iPawnW, iPawnH, wOff, hOff, ioVisit, ioMoves);
        }
        break;
      }
      else {
        wOff+= wInc;
        hOff+= hInc;
      }
    }
  }
}


void JumpinPlayerAI::SortSubBoards(BoardState *ioBoard, const int iDepth, const int iMode) {
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


int JumpinPlayerAI::GetIdxBestSubBoard(BoardState *ioBoard, const int iDepth, const int iMode) {
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


// A best= +1, B best= -1, Equal= 0
int JumpinPlayerAI::CompareBoardPair(const BoardState *iBoardA, const BoardState *iBoardB, const int iDepth, const int iMode) {
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
