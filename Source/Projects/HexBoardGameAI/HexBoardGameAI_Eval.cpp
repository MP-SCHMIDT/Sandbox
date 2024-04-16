#include "HexBoardGameAI.hpp"


// Standard lib

// Algo headers
#include "Math/Field.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


bool HexBoardGameAI::IsRedTurn(const int iDepth) {
  return ((idxTurn + iDepth) % turnPeriod) < streakRed;
}


void HexBoardGameAI::ComputeBoardScore(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] ComputeBoardScore on a null board\n");

  // Check for win state
  std::vector<std::vector<int>> Label= Field::AllocField2D(nW, nH, 0);
  for (int h= 0; h < nH; h++)
    if (ioBoard->Pawns[0][h] == -1)
      Label[0][h]= -1;
  for (int w= 0; w < nW; w++)
    if (ioBoard->Pawns[w][0] == +1)
      Label[w][0]= +1;

  bool wasUpdated= true;
  while (wasUpdated) {
    wasUpdated= false;
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (Label[w][h] == 0) {
          if (w - 1 >= 0 && Label[w - 1][h] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h]) Label[w][h]= Label[w - 1][h];
          if (w - 1 >= 0 && h - 1 >= 0 && Label[w - 1][h - 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 1]) Label[w][h]= Label[w - 1][h - 1];
          if (h - 1 >= 0 && Label[w][h - 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w][h - 1]) Label[w][h]= Label[w][h - 1];
          if (w + 1 < nW && Label[w + 1][h] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h]) Label[w][h]= Label[w + 1][h];
          if (w + 1 < nW && h + 1 < nH && Label[w + 1][h + 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h + 1]) Label[w][h]= Label[w + 1][h + 1];
          if (h + 1 < nH && Label[w][h + 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w][h + 1]) Label[w][h]= Label[w][h + 1];
          if (Label[w][h] != 0)
            wasUpdated= true;
        }
      }
    }
  }

  bool isWinRed= false;
  for (int w= 0; w < nW; w++)
    if (Label[w][nH - 1] == +1)
      isWinRed= true;
  bool isWinBlu= false;
  for (int h= 0; h < nH; h++)
    if (Label[nW - 1][h] == -1)
      isWinBlu= true;

  if (isWinRed) ioBoard->Score= +INT_MAX;
  else if (isWinBlu) ioBoard->Score= -INT_MAX;
  else {
    ioBoard->Score= 0;

    if (D.UI[ValEdgeConnect__].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0) {
            if (w > 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h]) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValEdgeConnect__].I();
            if (h > 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w][h - 1]) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValEdgeConnect__].I();
            if (w > 0 && h > 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 1]) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValEdgeConnect__].I();
          }
        }
      }
    }

    if (D.UI[ValCornConnect__].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0) {
            if (w > 1 && h > 0 &&
                ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 2][h - 1] &&
                ioBoard->Pawns[w - 1][h] == 0 &&
                ioBoard->Pawns[w - 1][h - 1] == 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValCornConnect__].I();
            if (w > 0 && h > 1 &&
                ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 2] &&
                ioBoard->Pawns[w][h - 1] == 0 &&
                ioBoard->Pawns[w - 1][h - 1] == 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValCornConnect__].I();
            if (w < nW - 1 && h > 0 &&
                ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h - 1] &&
                ioBoard->Pawns[w + 1][h] == 0 &&
                ioBoard->Pawns[w][h - 1] == 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValCornConnect__].I();
          }
        }
      }
    }
  }
}
