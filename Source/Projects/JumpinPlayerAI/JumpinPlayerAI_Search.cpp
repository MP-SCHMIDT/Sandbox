#include "JumpinPlayerAI.hpp"


// Standard lib

// Algo headers
#include "Math/Field.hpp"
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
  for (BoardState *subBoard : RootBoard->SubBoards)
    DeleteBoard(subBoard);
  RootBoard->SubBoards.clear();
  RootBoard->Move= {-1, -1, -1, -1};
  ComputeBoardScore(RootBoard);

  // TODO add alpha beta pruning

  // Run the recursive search with iterative deepening
  for (int iterMaxSearchDepth= 1; iterMaxSearchDepth <= D.UI[MaxSearchDepth__].I(); iterMaxSearchDepth++)
    RecursiveTreeSearch(RootBoard, 0, iterMaxSearchDepth);

  // Stop the timer
  thinkTime= Timer::PopTimer();
}


void JumpinPlayerAI::RecursiveTreeSearch(BoardState *ioBoard, const int iDepth, const int iMaxDepth) {
  if (ioBoard == nullptr) printf("[ERROR] RecursiveTreeSearch on a null board\n");

  // Handle leaf board
  if (iDepth >= iMaxDepth ||
      ioBoard->Score == +INT_MAX ||
      ioBoard->Score == -INT_MAX ||
      (D.UI[MaxThinkTime____].D() > 0.0 && Timer::CheckTimer() >= D.UI[MaxThinkTime____].D()) ||
      (D.UI[MaxTreeBoards___].I() > 0 && nbTreeBoards >= D.UI[MaxTreeBoards___].I())) {
    ioBoard->NashScore= ioBoard->Score;
    ioBoard->NashNbSteps= 0;
    return;
  }

  // Create and score the sub boards if needed
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
    // Add the turn skip move if no other valid move
    if (Moves.empty())
      Moves.push_back(std::array<int, 4>{-1, -1, -1, -1});
    // Create the sub board for each move
    for (std::array<int, 4> move : Moves) {
      BoardState *newBoard= CreateBoard(ioBoard->Pawns, move, iDepth + 1);
      if (move[0] != -1 && move[1] != -1 && move[2] != -1 && move[3] != -1) {
        newBoard->Pawns[move[2]][move[3]]= newBoard->Pawns[move[0]][move[1]];
        newBoard->Pawns[move[0]][move[1]]= 0;
      }
      ioBoard->SubBoards.push_back(newBoard);
      nbTreeBoards++;
    }
  }

  // Recursively search the sub boards and update Nash
  for (BoardState *subBoard : ioBoard->SubBoards) {
    ComputeBoardScore(subBoard);
    RecursiveTreeSearch(subBoard, iDepth + 1, iMaxDepth);
  }
  SortSubBoards(ioBoard, iDepth);
  ioBoard->NashScore= ioBoard->SubBoards[0]->NashScore;
  ioBoard->NashNbSteps= ioBoard->SubBoards[0]->NashNbSteps + 1;

  return;
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


void JumpinPlayerAI::SortSubBoards(BoardState *ioBoard, const int iDepth) {
  if (ioBoard == nullptr) printf("[ERROR] SortSubBoards on a null board\n");

  // Bubble sort of the sub boards in Nash order
  for (int k0= 0; k0 < (int)ioBoard->SubBoards.size(); k0++) {
    for (int k1= k0 + 1; k1 < (int)ioBoard->SubBoards.size(); k1++) {
      if ((IsRedTurn(iDepth) && ioBoard->SubBoards[k0]->NashScore < ioBoard->SubBoards[k1]->NashScore) ||
          (!IsRedTurn(iDepth) && ioBoard->SubBoards[k0]->NashScore > ioBoard->SubBoards[k1]->NashScore) ||
          (ioBoard->SubBoards[k0]->NashScore == ioBoard->SubBoards[k1]->NashScore &&
           ioBoard->SubBoards[k0]->NashNbSteps > ioBoard->SubBoards[k1]->NashNbSteps)) {
        BoardState *tmp= ioBoard->SubBoards[k0];
        ioBoard->SubBoards[k0]= ioBoard->SubBoards[k1];
        ioBoard->SubBoards[k1]= tmp;
      }
    }
  }
}
