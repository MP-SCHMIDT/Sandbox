#include "JumpinPlayerAI.hpp"


// Standard lib

// Algo headers

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


bool JumpinPlayerAI::IsRedTurn(const int iDepth) {
  return (idxTurn + iDepth) % 2 == 0;
}


void JumpinPlayerAI::ComputeBoardScore(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] ComputeBoardScore on a null board\n");

  // Check for win state
  bool isWinRed= true;
  bool isWinBlu= true;
  for (int w= 0; w < nW && isWinRed; w++)
    for (int h= nH - D.UI[StartingRows____].I(); h < nH && isWinRed; h++)
      if (ioBoard->Pawns[w][h] <= 0) isWinRed= false;
  for (int w= 0; w < nW && isWinBlu; w++)
    for (int h= 0; h < D.UI[StartingRows____].I() && isWinBlu; h++)
      if (ioBoard->Pawns[w][h] >= 0) isWinBlu= false;

  if (isWinRed) ioBoard->Score= +INT_MAX;
  else if (isWinBlu) ioBoard->Score= -INT_MAX;
  else {
    // Reset the score
    ioBoard->Score= 0;

    // Add score for total pawn advance
    if (D.UI[ValPushTotal____].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] > 0) ioBoard->Score+= (h + 1) * D.UI[ValPushTotal____].I();
          if (ioBoard->Pawns[w][h] < 0) ioBoard->Score-= (nH - h) * D.UI[ValPushTotal____].I();
        }
      }
    }

    // Add score for lagging pawn advance
    if (D.UI[ValPushLast_____].I() != 0) {
      bool foundLastRed= false;
      for (int h= 0; h < nH && !foundLastRed; h++) {
        for (int w= 0; w < nW && !foundLastRed; w++) {
          if (ioBoard->Pawns[w][h] > 0) {
            ioBoard->Score+= (h + 1) * D.UI[ValPushLast_____].I();
            foundLastRed= true;
          }
        }
      }
      bool foundLastBlu= false;
      for (int h= nH - 1; h >= 0 && !foundLastBlu; h--) {
        for (int w= 0; w < nW && !foundLastBlu; w++) {
          if (ioBoard->Pawns[w][h] < 0) {
            ioBoard->Score-= (nH - h) * D.UI[ValPushLast_____].I();
            foundLastBlu= true;
          }
        }
      }
    }

    // Add score for soft stranded pawns
    if (D.UI[ValSoftStranded_].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0 &&
              (w - 1 < 0 || ioBoard->Pawns[w - 1][h] != ioBoard->Pawns[w][h]) &&
              (w + 1 >= nW || ioBoard->Pawns[w + 1][h] != ioBoard->Pawns[w][h]) &&
              (h - 1 < 0 || ioBoard->Pawns[w][h - 1] != ioBoard->Pawns[w][h]) &&
              (h + 1 >= nH || ioBoard->Pawns[w][h + 1] != ioBoard->Pawns[w][h]))
            ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValSoftStranded_].I();
        }
      }
    }

    // Add score for hard stranded pawns
    if (D.UI[ValHardStranded_].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0 &&
              (w - 1 < 0 || ioBoard->Pawns[w - 1][h] == 0) &&
              (w + 1 >= nW || ioBoard->Pawns[w + 1][h] == 0) &&
              (h - 1 < 0 || ioBoard->Pawns[w][h - 1] == 0) &&
              (h + 1 >= nH || ioBoard->Pawns[w][h + 1] == 0))
            ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ValHardStranded_].I();
        }
      }
    }
  }
}
