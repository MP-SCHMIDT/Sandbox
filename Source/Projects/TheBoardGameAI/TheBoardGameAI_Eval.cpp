#include "TheBoardGameAI.hpp"


// Standard lib

// Algo headers
#include "Math/Field.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


bool TheBoardGameAI::IsRedTurn(const int iDepth) {
  return ((idxTurn + iDepth) % (streakRed + streakBlu)) < streakRed;
}


void TheBoardGameAI::ComputeBoardScoreHex(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] ComputeBoardScoreHex on a null board\n");

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
          if (w - 1 >= 0 && Label[w - 1][h] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h])
            Label[w][h]= Label[w - 1][h];
          if (w - 1 >= 0 && h - 1 >= 0 && Label[w - 1][h - 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 1])
            Label[w][h]= Label[w - 1][h - 1];
          if (h - 1 >= 0 && Label[w][h - 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w][h - 1])
            Label[w][h]= Label[w][h - 1];
          if (w + 1 < nW && Label[w + 1][h] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h])
            Label[w][h]= Label[w + 1][h];
          if (w + 1 < nW && h + 1 < nH && Label[w + 1][h + 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h + 1])
            Label[w][h]= Label[w + 1][h + 1];
          if (h + 1 < nH && Label[w][h + 1] != 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w][h + 1])
            Label[w][h]= Label[w][h + 1];
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

    if (D.UI[HexEdgeConnect__].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0) {
            if (w > 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h])
              ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[HexEdgeConnect__].I();
            if (h > 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w][h - 1])
              ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[HexEdgeConnect__].I();
            if (w > 0 && h > 0 && ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 1])
              ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[HexEdgeConnect__].I();
          }
        }
      }
    }

    if (D.UI[HexCornConnect__].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0) {
            if (w > 1 && h > 0 &&
                ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 2][h - 1] &&
                ioBoard->Pawns[w - 1][h] == 0 &&
                ioBoard->Pawns[w - 1][h - 1] == 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[HexCornConnect__].I();
            if (w > 0 && h > 1 &&
                ioBoard->Pawns[w][h] == ioBoard->Pawns[w - 1][h - 2] &&
                ioBoard->Pawns[w][h - 1] == 0 &&
                ioBoard->Pawns[w - 1][h - 1] == 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[HexCornConnect__].I();
            if (w < nW - 1 && h > 0 &&
                ioBoard->Pawns[w][h] == ioBoard->Pawns[w + 1][h - 1] &&
                ioBoard->Pawns[w + 1][h] == 0 &&
                ioBoard->Pawns[w][h - 1] == 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[HexCornConnect__].I();
          }
        }
      }
    }
  }
}


void TheBoardGameAI::ComputeBoardScoreJmp(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] ComputeBoardScoreJmp on a null board\n");

  // Check for win state
  bool isWinRed= true;
  for (int w= 0; w < nW && isWinRed; w++)
    for (int h= std::max(0, nH - 2); h < nH && isWinRed; h++)
      if (ioBoard->Pawns[w][h] <= 0) isWinRed= false;
  bool isWinBlu= true;
  for (int w= 0; w < nW && isWinBlu; w++)
    for (int h= 0; h < std::min(2, nH) && isWinBlu; h++)
      if (ioBoard->Pawns[w][h] >= 0) isWinBlu= false;

  if (isWinRed) ioBoard->Score= +INT_MAX;
  else if (isWinBlu) ioBoard->Score= -INT_MAX;
  else {
    // Reset the score
    ioBoard->Score= 0;

    // Add score for total pawn advance
    if (D.UI[JmpPushTotal____].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] > 0) ioBoard->Score+= (h + 1) * D.UI[JmpPushTotal____].I();
          if (ioBoard->Pawns[w][h] < 0) ioBoard->Score-= (nH - h) * D.UI[JmpPushTotal____].I();
        }
      }
    }

    // Add score for lagging pawn advance
    // TODO change for sum of penal on back row
    if (D.UI[JmpPushLast_____].I() != 0) {
      bool foundLastRed= false;
      for (int h= 0; h < nH && !foundLastRed; h++) {
        for (int w= 0; w < nW && !foundLastRed; w++) {
          if (ioBoard->Pawns[w][h] > 0) {
            ioBoard->Score+= (h + 1) * D.UI[JmpPushLast_____].I();
            foundLastRed= true;
          }
        }
      }
      bool foundLastBlu= false;
      for (int h= nH - 1; h >= 0 && !foundLastBlu; h--) {
        for (int w= 0; w < nW && !foundLastBlu; w++) {
          if (ioBoard->Pawns[w][h] < 0) {
            ioBoard->Score-= (nH - h) * D.UI[JmpPushLast_____].I();
            foundLastBlu= true;
          }
        }
      }
    }

    // Add score for soft stranded pawns
    if (D.UI[JmpSoftStranded_].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0 &&
              (w - 1 < 0 || ioBoard->Pawns[w - 1][h] != ioBoard->Pawns[w][h]) &&
              (w + 1 >= nW || ioBoard->Pawns[w + 1][h] != ioBoard->Pawns[w][h]) &&
              (h - 1 < 0 || ioBoard->Pawns[w][h - 1] != ioBoard->Pawns[w][h]) &&
              (h + 1 >= nH || ioBoard->Pawns[w][h + 1] != ioBoard->Pawns[w][h]))
            ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[JmpSoftStranded_].I();
        }
      }
    }

    // Add score for hard stranded pawns
    if (D.UI[JmpHardStranded_].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] != 0 &&
              (w - 1 < 0 || ioBoard->Pawns[w - 1][h] == 0) &&
              (w + 1 >= nW || ioBoard->Pawns[w + 1][h] == 0) &&
              (h - 1 < 0 || ioBoard->Pawns[w][h - 1] == 0) &&
              (h + 1 >= nH || ioBoard->Pawns[w][h + 1] == 0))
            ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[JmpHardStranded_].I();
        }
      }
    }
  }
}

void TheBoardGameAI::ComputeBoardScoreChk(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] ComputeBoardScoreChk on a null board\n");

  // Check for win state
  bool isWinRed= true;
  bool isWinBlu= true;
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (ioBoard->Pawns[w][h] < 0) isWinRed= false;
      if (ioBoard->Pawns[w][h] > 0) isWinBlu= false;
    }
  }

  if (isWinRed) ioBoard->Score= +INT_MAX;
  else if (isWinBlu) ioBoard->Score= -INT_MAX;
  else {
    // Reset the score
    ioBoard->Score= 0;

    // Count the material
    for (int w= 0; w < nW; w++)
      for (int h= 0; h < nH; h++)
        if (ioBoard->Pawns[w][h] != 0) ioBoard->Score+= ioBoard->Pawns[w][h] * D.UI[ChkMaterial_____].I();
  }
}