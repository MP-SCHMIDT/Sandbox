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
        ioMoves.push_back(std::vector<std::array<int, 2>>{{w, h}});
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
        for (std::array<int, 2> destination : Destinations)
          ioMoves.push_back(std::vector<std::array<int, 2>>{{w, h}, destination});
      }
    }
  }
}


void TheBoardGameAI::FindPossibleMovesChk(BoardState *ioBoard, const int iDepth, std::vector<std::vector<std::array<int, 2>>> &ioMoves) {
  bool hasFoundCaptures= false;
  // Sweep through the current player pawns
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if ((ioBoard->Pawns[w][h] > 0 && IsRedTurn(iDepth)) || (ioBoard->Pawns[w][h] < 0 && !IsRedTurn(iDepth))) {
        // Compute possible captures for the current pawn
        for (int idxDir= 0; idxDir < 4; idxDir++) {
          const int wHalf= (idxDir == 0 || idxDir == 1) ? (w - 1) : (w + 1);
          const int hHalf= (idxDir == 1 || idxDir == 3) ? (h - 1) : (h + 1);
          const int wFull= (idxDir == 0 || idxDir == 1) ? (w - 2) : (w + 2);
          const int hFull= (idxDir == 1 || idxDir == 3) ? (h - 2) : (h + 2);
          if (wFull >= 0 && wFull < nW && hFull >= 0 && hFull < nH &&
              ioBoard->Pawns[wFull][hFull] == 0 && ioBoard->Pawns[wHalf][hHalf] != 0 &&
              ioBoard->Pawns[w][h] != ioBoard->Pawns[wHalf][hHalf]) {
            hasFoundCaptures= true;
            ioMoves.push_back(std::vector<std::array<int, 2>>{{w, h}, {wHalf, hHalf}, {wFull, hFull}});
          }
        }
      }
    }
  }
  // Sweep through the current player pawns
  if (!hasFoundCaptures) {
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if ((ioBoard->Pawns[w][h] > 0 && IsRedTurn(iDepth)) || (ioBoard->Pawns[w][h] < 0 && !IsRedTurn(iDepth))) {
          // Compute possible moves for the current pawn
          for (int idxDir= 0; idxDir < 4; idxDir++) {
            const int wHalf= (idxDir == 0 || idxDir == 1) ? (w - 1) : (w + 1);
            const int hHalf= (idxDir == 1 || idxDir == 3) ? (h - 1) : (h + 1);
            if (wHalf >= 0 && wHalf < nW && hHalf >= 0 && hHalf < nH &&
                ioBoard->Pawns[wHalf][hHalf] == 0) {
              ioMoves.push_back(std::vector<std::array<int, 2>>{{w, h}, {wHalf, hHalf}});
            }
          }
        }
      }
    }
  }

  // TODO implement move with recursive jumps and captures
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
          ioDestinations.push_back({wOff, hOff});
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


void TheBoardGameAI::ExecuteMoveHex(BoardState *ioBoard, const int iDepth) {
  ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]]= IsRedTurn(iDepth) ? +1 : -1;
}


void TheBoardGameAI::ExecuteMoveJmp(BoardState *ioBoard, const int iDepth) {
  (void)iDepth;
  ioBoard->Pawns[ioBoard->Move[1][0]][ioBoard->Move[1][1]]= ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]];
  ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]]= 0;
}


void TheBoardGameAI::ExecuteMoveChk(BoardState *ioBoard, const int iDepth) {
  (void)iDepth;
  if (ioBoard->Move.size() == 3) {
    ioBoard->Pawns[ioBoard->Move[2][0]][ioBoard->Move[2][1]]= ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]];
    ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]]= 0;
    ioBoard->Pawns[ioBoard->Move[1][0]][ioBoard->Move[1][1]]= 0;
  }
  else {
    ioBoard->Pawns[ioBoard->Move[1][0]][ioBoard->Move[1][1]]= ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]];
    ioBoard->Pawns[ioBoard->Move[0][0]][ioBoard->Move[0][1]]= 0;
  }
}
