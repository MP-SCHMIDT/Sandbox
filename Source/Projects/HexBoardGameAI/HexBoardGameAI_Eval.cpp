#include "HexBoardGameAI.hpp"


// Standard lib

// Algo headers
#include "Math/Field.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


bool HexBoardGameAI::IsRedTurn(const int iDepth) {
  // TODO test turn swap imbalence to see (e.g. using %3).
  // could help with tournament testing by giving edge ton one player and limiting possibility of antigame strat
  return (idxTurn + iDepth) % 2 == 0;
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
          if (w - 1 >= 0 &&               Label[w - 1][h    ] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h    ]) { Label[w][h]= Label[w - 1][h    ]; wasUpdated= true; }
          if (w - 1 >= 0 && h - 1 >= 0 && Label[w - 1][h - 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 1]) { Label[w][h]= Label[w - 1][h - 1]; wasUpdated= true; }
          if (              h - 1 >= 0 && Label[w    ][h - 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w    ][h - 1]) { Label[w][h]= Label[w    ][h - 1]; wasUpdated= true; }
          if (w + 1 < nW &&               Label[w + 1][h    ] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h    ]) { Label[w][h]= Label[w + 1][h    ]; wasUpdated= true; }
          if (w + 1 < nW && h + 1 < nH && Label[w + 1][h + 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h + 1]) { Label[w][h]= Label[w + 1][h + 1]; wasUpdated= true; }
          if (              h + 1 < nH && Label[w    ][h + 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w    ][h + 1]) { Label[w][h]= Label[w    ][h + 1]; wasUpdated= true; }
        }
      }
    }
  }

  bool isWinRed= false;
  for (int w= 0; w < nW; w++)
    if (Label[w][nH-1] == +1)
      isWinRed= true;
  bool isWinBlu= false;
  for (int h= 0; h < nH; h++)
    if (Label[nW-1][h] == -1)
      isWinBlu= true;

  if (isWinRed) ioBoard->Score= +INT_MAX;
  else if (isWinBlu) ioBoard->Score= -INT_MAX;
  else ioBoard->Score= 0;
}
