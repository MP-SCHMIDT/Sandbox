#include "TheBoardGameAI.hpp"


// Standard lib

// Algo headers
#include "Math/Field.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void TheBoardGameAI::FindPossibleMovesHex(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves) {
  (void)iDepth;
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (ioBoard->Pawns[w][h] == 0) {
        ioMoves.push_back(std::vector<std::array<int, 2>>(1, std::array<int, 2>{w, h}));
      }
    }
  }
}


void TheBoardGameAI::FindPossibleMovesJmp(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves) {
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if ((ioBoard->Pawns[w][h] > 0 && IsRedTurn(iDepth)) || (ioBoard->Pawns[w][h] < 0 && !IsRedTurn(iDepth))) {
        // Compute possible moves for the current pawn
        std::vector<std::vector<bool>> Visit= Field::AllocField2D(nW, nH, false);
        std::vector<std::array<int, 2>> Destinations;
        Visit[w][h]= true;
        const int oldPawn= ioBoard->Pawns[w][h];
        ioBoard->Pawns[w][h]= 0;
        JmpRecursivePawnMoves(ioBoard, w, h, Visit, Destinations);
        ioBoard->Pawns[w][h]= oldPawn;
        if (!Destinations.empty()) {
          ioMoves.push_back(std::vector<std::array<int, 2>>(1, std::array<int, 2>{w, h}));
          for (std::array<int, 2> destination : Destinations)
            ioMoves[ioMoves.size() - 1].push_back(destination);
        }
      }
    }
  }
  if (ioMoves.empty())
    ioMoves.push_back(std::vector<std::array<int, 2>>(2, std::array<int, 2>{-1, -1}));
}


void TheBoardGameAI::FindPossibleMovesChk(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves) {
  (void)ioBoard;
  (void)iDepth;
  (void)ioMoves;

  // TODO implement move with recursive jumps and captues
}


void TheBoardGameAI::JmpRecursivePawnMoves(BoardState *ioBoard,
                                           const int iJumpW, const int iJumpH,
                                           std::vector<std::vector<bool>> &ioVisit,
                                           std::vector<std::array<int, 2>> &ioDestinations) {
  if (ioBoard == nullptr) printf("[ERROR] JmpRecursivePawnMoves on a null board\n");

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
          ioDestinations.push_back(std::array<int, 2>{wOff, hOff});
          JmpRecursivePawnMoves(ioBoard, wOff, hOff, ioVisit, ioDestinations);
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
