#include "JumpinPlayerAI.hpp"


// Standard lib
#include <cmath>
#include <numbers>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/Bresenham.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
JumpinPlayerAI::JumpinPlayerAI() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void JumpinPlayerAI::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("BoardW__________", 4));
    D.UI.push_back(ParamUI("BoardH__________", 5));
    D.UI.push_back(ParamUI("StartingRows____", 2));
    D.UI.push_back(ParamUI("SinglePlayer____", 1));
    D.UI.push_back(ParamUI("SearchDepth_____", 2));
    D.UI.push_back(ParamUI("TreeStepDist____", 0.3));
    D.UI.push_back(ParamUI("TreeStepRadians_", 0.5));
    D.UI.push_back(ParamUI("TreeFactDist____", 0.5));
    D.UI.push_back(ParamUI("TreeFactRadians_", 0.125));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.e-2));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));
  }

  if (D.UI.size() != VerboseLevel____ + 1) {
    printf("[ERROR] Invalid parameter count in UI\n");
  }

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool JumpinPlayerAI::CheckAlloc() {
  if (D.UI[BoardW__________].hasChanged()) isAllocated= false;
  if (D.UI[BoardH__________].hasChanged()) isAllocated= false;
  if (D.UI[StartingRows____].hasChanged()) isAllocated= false;
  if (D.UI[SinglePlayer____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool JumpinPlayerAI::CheckRefresh() {
  if (D.UI[SearchDepth_____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void JumpinPlayerAI::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Clear any existing board tree
  if (RootBoard != nullptr) DeleteBoard(RootBoard);

  // Get UI parameters
  nW= std::max(D.UI[BoardW__________].I(), 1);
  nH= std::max(D.UI[BoardH__________].I(), 1);
  idxTurn= 0;
  wSel= -1;
  hSel= -1;

  // Set box dimensions
  const double voxSize= 1.0 / (double)nH;
  D.boxMin= {0.5 - 0.5 * voxSize,
             0.5 - 0.5 * voxSize * double(nW),
             0.5 - 0.5 * voxSize * double(nH)};
  D.boxMax= {0.5 + 0.5 * voxSize,
             0.5 + 0.5 * voxSize * double(nW),
             0.5 + 0.5 * voxSize * double(nH)};

  // Create and initialize root board with pawns
  RootBoard= CreateBoard();
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (h < D.UI[StartingRows____].I()) RootBoard->Pawns[w][h]= +1;
      if (!D.UI[SinglePlayer____].B())
        if (h >= nH - D.UI[StartingRows____].I()) RootBoard->Pawns[w][h]= -1;
    }
  }
}


// Refresh the project
void JumpinPlayerAI::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Compute root board score and moves
  ComputeBoardScore(RootBoard);
  ComputeBoardMoves(RootBoard, 0);
  ConsolidateBoardScores(RootBoard, 0);

  if (D.Plot.size() < 1) D.Plot.resize(1);
  D.Plot[0].name= "Score";
  D.Plot[0].isSymmetric= true;
  D.Plot[0].val.push_back(RootBoard->Score);
}


// Handle keypress
void JumpinPlayerAI::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Get cursor coordinates on the board
  const int wCursor= (int)std::floor((double)nW * (D.mouseProjX[1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
  const int hCursor= (int)std::floor((double)nH * (D.mouseProjX[2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));

  // Handle keyboard action
  if (key == 'E' || key == 'W' || key == 'B') {
    if (wCursor >= 0 && wCursor < nW && hCursor >= 0 && hCursor < nH) {
      if (key == 'E') RootBoard->Pawns[wCursor][hCursor]= 0;
      if (key == 'W') RootBoard->Pawns[wCursor][hCursor]= +1;
      if (key == 'B') RootBoard->Pawns[wCursor][hCursor]= -1;
      ComputeBoardScore(RootBoard);
      ComputeBoardMoves(RootBoard, 0);
      ConsolidateBoardScores(RootBoard, 0);
    }
  }

  if (key == 'S') {
    wSel= wCursor;
    hSel= hCursor;
  }

  if (key == 'D') {
    std::array<int, 2> target= {wCursor, hCursor};
    if (wSel >= 0 && wSel < nW && hSel >= 0 && hSel < nH) {
      for (std::array<int, 2> move : RootBoard->Moves[wSel][hSel]) {
        if (move == target) {
          RootBoard->Pawns[target[0]][target[1]]= RootBoard->Pawns[wSel][hSel];
          RootBoard->Pawns[wSel][hSel]= 0;
          break;
        }
      }
    }
    ComputeBoardScore(RootBoard);
    ComputeBoardMoves(RootBoard, 0);
    ConsolidateBoardScores(RootBoard, 0);
    idxTurn++;
  }

  if (key == 'A' || key == 'Z') {
    bool isSet= false;
    int scoreBest= 0, wBest= 0, hBest= 0, kBest= 0;
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        for (int k= 0; k < (int)RootBoard->Moves[w][h].size(); k++) {
          bool newBest= false;
          if (key == 'A' && RootBoard->Pawns[w][h] > 0 && (!isSet || scoreBest < RootBoard->SubBoards[w][h][k]->Score)) newBest= true;
          if (key == 'Z' && RootBoard->Pawns[w][h] < 0 && (!isSet || scoreBest > RootBoard->SubBoards[w][h][k]->Score)) newBest= true;
          if (newBest) {
            isSet= true;
            scoreBest= RootBoard->SubBoards[w][h][k]->Score;
            wBest= w;
            hBest= h;
            kBest= k;
          }
        }
      }
    }
    if (isSet) {
      std::array<int, 2> moveBest= RootBoard->Moves[wBest][hBest][kBest];
      RootBoard->Pawns[moveBest[0]][moveBest[1]]= RootBoard->Pawns[wBest][hBest];
      RootBoard->Pawns[wBest][hBest]= 0;
    }
    ComputeBoardScore(RootBoard);
    ComputeBoardMoves(RootBoard, 0);
    ConsolidateBoardScores(RootBoard, 0);
    idxTurn++;
  }

  if (key == 'F' || key == 'G') {
    std::array<int, 2> beg, end;
    if (key == 'F') {
      beg= RootBoard->MoveBestWhite[0];
      end= RootBoard->MoveBestWhite[1];
      RootBoard->MoveBestWhite= {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, -1}};
    }
    else {
      beg= RootBoard->MoveBestBlack[0];
      end= RootBoard->MoveBestBlack[1];
      RootBoard->MoveBestBlack= {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, -1}};
    }
    if (beg[0] != -1 && beg[1] != -1 && end[0] != -1 && end[1] != -1) {
      RootBoard->Pawns[end[0]][end[1]]= RootBoard->Pawns[beg[0]][beg[1]];
      RootBoard->Pawns[beg[0]][beg[1]]= 0;
      idxTurn++;
      ComputeBoardScore(RootBoard);
      ComputeBoardMoves(RootBoard, 0);
      ConsolidateBoardScores(RootBoard, 0);
    }
  }

  if (D.Plot.size() < 1) D.Plot.resize(1);
  D.Plot[0].name= "Score";
  D.Plot[0].isSymmetric= true;
  D.Plot[0].val.push_back(RootBoard->Score);
}


// Animate the project
void JumpinPlayerAI::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
}


// Draw the project
void JumpinPlayerAI::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  const float voxSize= 1.0 / (float)nH;

  // Draw the board
  if (D.displayMode1) {
    glPushMatrix();
    glTranslatef(0.5f, 0.0f, 0.0f);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if ((w + h) % 2 == 0) glColor3f(0.45f, 0.45f, 0.45f);
        else glColor3f(0.55f, 0.55f, 0.55f);
        glRectf(D.boxMin[1] + float(w) * voxSize, D.boxMin[2] + float(h) * voxSize,
                D.boxMin[1] + float(w + 1) * voxSize, D.boxMin[2] + float(h + 1) * voxSize);
      }
    }
    glPopMatrix();
  }

  // Draw the pawns
  if (D.displayMode2) {
    glEnable(GL_LIGHTING);
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (RootBoard->Pawns[w][h] != 0) {
          float r= 0.2f, g= 0.2f, b= 0.2f;
          if (RootBoard->Pawns[w][h] < 0) r= g= b= 0.2f;
          if (RootBoard->Pawns[w][h] > 0) r= g= b= 0.8f;
          if (w == wSel && h == hSel) g+= 0.2f;
          glColor3f(r, g, b);
          glPushMatrix();
          glTranslatef(D.boxMin[0] + 0.5 * voxSize,
                       D.boxMin[1] + 0.5 * voxSize + float(w) * voxSize,
                       D.boxMin[2] + 0.5 * voxSize + float(h) * voxSize);
          glutSolidSphere(0.45 * voxSize, 36, 10);
          glPopMatrix();
        }
      }
    }
    glDisable(GL_LIGHTING);
  }

  // Draw the available moves for the selected pawn
  if (D.displayMode3) {
    glLineWidth(2.0f);
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    if (wSel >= 0 && wSel < nW && hSel >= 0 && hSel < nH) {
      for (int k= 0; k < (int)RootBoard->Moves[wSel][hSel].size(); k++) {
        std::array<int, 2> move= RootBoard->Moves[wSel][hSel][k];
        BoardState *subBoard= RootBoard->SubBoards[wSel][hSel][k];
        float r, g, b;
        const float relScore= (float)(subBoard->Score - RootBoard->Score);
        Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * relScore, r, g, b);
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(0.0f, (float)move[0], (float)move[1]);
        glutWireCube(0.95f);
        glPopMatrix();
      }
    }
    glPopMatrix();
    glLineWidth(1.0f);
  }

  // Draw the board tree
  if (D.displayMode4) {
    float px= 0.5f * (D.boxMin[0] + D.boxMax[0]);
    float py= D.boxMax[1] + 0.6f * (D.boxMax[1] - D.boxMin[1]);
    float pz= D.boxMin[2] + 0.5f * (D.boxMax[2] - D.boxMin[2]);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    DrawBoardTree(RootBoard, 1, px, py, pz, 0.0f, 0.0f,
                  D.UI[TreeStepDist____].F(), D.UI[TreeStepRadians_].F(),
                  D.UI[TreeFactDist____].F(), D.UI[TreeFactRadians_].F());
    glEnd();
    glLineWidth(1.0f);
    glPointSize(10.0f);
    glBegin(GL_POINTS);
    DrawBoardTree(RootBoard, 0, px, py, pz, 0.0f, 0.0f,
                  D.UI[TreeStepDist____].F(), D.UI[TreeStepRadians_].F(),
                  D.UI[TreeFactDist____].F(), D.UI[TreeFactRadians_].F());
    glEnd();
    glPointSize(1.0f);
  }
}


void JumpinPlayerAI::DrawBoardTree(const BoardState *iBoard, const int iDrawMode,
                                   const float px, const float py, const float pz,
                                   const float dist, const float radians,
                                   const float distStep, const float radiansStep,
                                   const float distFact, const float radiansFact) {
  if (iBoard == nullptr) printf("Error: Drawing a null board");

  float distOff= dist + distStep;
  float radiansOff= radians;
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      for (int k= 0; k < (int)iBoard->Moves[w][h].size(); k++) {
        BoardState *subBoard= iBoard->SubBoards[w][h][k];
        if (iDrawMode == 0) {
          float r, g, b;
          const float relScore= (float)(subBoard->ScoreBestWhite - RootBoard->Score);
          Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * relScore, r, g, b);
          glColor3f(r, g, b);
          glVertex3f(px + 0.01f, py + distOff * std::cos(radiansOff), pz + distOff * std::sin(radiansOff));
        }
        else {
          if (iBoard->MoveBestWhite[0] == std::array<int, 2>{w, h} && iBoard->MoveBestWhite[1] == iBoard->Moves[w][h][k]) {
            glColor3f(1.0f, 1.0f, 1.0f);
            glVertex3f(px + 0.01f, py + dist * std::cos(radians), pz + dist * std::sin(radians));
            glVertex3f(px + 0.01f, py + distOff * std::cos(radiansOff), pz + distOff * std::sin(radiansOff));
          }
          if (iBoard->MoveBestBlack[0] == std::array<int, 2>{w, h} && iBoard->MoveBestBlack[1] == iBoard->Moves[w][h][k]) {
            glColor3f(0.0f, 0.0f, 0.0f);
            glVertex3f(px + 0.01f, py + dist * std::cos(radians), pz + dist * std::sin(radians));
            glVertex3f(px + 0.01f, py + distOff * std::cos(radiansOff), pz + distOff * std::sin(radiansOff));
          }
          float r, g, b;
          const float relScore= (float)(subBoard->Score - RootBoard->Score);
          Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * relScore, r, g, b);
          glColor3f(r, g, b);
          glVertex3f(px, py + dist * std::cos(radians), pz + dist * std::sin(radians));
          glVertex3f(px, py + distOff * std::cos(radiansOff), pz + distOff * std::sin(radiansOff));
        }
        DrawBoardTree(subBoard, iDrawMode, px, py, pz, distOff, radiansOff, distStep * distFact, radiansStep * radiansFact, distFact, radiansFact);
        radiansOff+= radiansStep;
      }
    }
  }
}


JumpinPlayerAI::BoardState *JumpinPlayerAI::CreateBoard() {
  BoardState *newBoard= new BoardState;
  newBoard->Pawns= Field::AllocField2D(nW, nH, 0);
  newBoard->Moves= Field::AllocField2D(nW, nH, std::vector<std::array<int, 2>>());
  newBoard->SubBoards= Field::AllocField2D(nW, nH, std::vector<BoardState *>());
  newBoard->Score= 0;
  newBoard->ScoreBestWhite= 0;
  newBoard->ScoreBestBlack= 0;
  newBoard->StepsBestWhite= 0;
  newBoard->StepsBestBlack= 0;
  newBoard->MoveBestWhite= {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, -1}};
  newBoard->MoveBestBlack= {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, -1}};
  return newBoard;
}


void JumpinPlayerAI::DeleteBoard(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("Error: Deleting a null board");

  for (int w= 0; w < nW; w++)
    for (int h= 0; h < nH; h++)
      for (BoardState *subBoard : ioBoard->SubBoards[w][h])
        DeleteBoard(subBoard);
  delete ioBoard;
  ioBoard= nullptr;
}


void JumpinPlayerAI::DeleteSubBoards(BoardState *ioBoard, const int w, const int h) {
  if (ioBoard == nullptr) printf("Error: Deleting subBoards in a null board");

  for (BoardState *subBoard : ioBoard->SubBoards[w][h]) {
    DeleteBoard(subBoard);
  }
  ioBoard->SubBoards[w][h].clear();
}


void JumpinPlayerAI::ComputeBoardScore(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("Error: Computing score of a null board");

  // Reset score
  ioBoard->Score= 0;

  // Add score for pawn advance
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (ioBoard->Pawns[w][h] > 0) ioBoard->Score+= (h + 1) * 10;
      if (ioBoard->Pawns[w][h] < 0) ioBoard->Score-= (nH - h) * 10;
    }
  }
}


void JumpinPlayerAI::ConsolidateBoardScores(BoardState *ioBoard, const int iDepth) {
  if (ioBoard == nullptr) printf("Error: Consolidating score of a null board");

  ioBoard->ScoreBestWhite= ioBoard->Score;
  ioBoard->ScoreBestBlack= ioBoard->Score;
  ioBoard->StepsBestWhite= 0;
  ioBoard->StepsBestBlack= 0;
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      for (int k= 0; k < (int)ioBoard->Moves[w][h].size(); k++) {
        BoardState *subBoard= ioBoard->SubBoards[w][h][k];
        ConsolidateBoardScores(subBoard, iDepth + 1);
        if (D.UI[SinglePlayer____].B()) {
          if (ioBoard->Pawns[w][h] > 0) {
            if (ioBoard->ScoreBestWhite < subBoard->ScoreBestWhite ||
                (ioBoard->ScoreBestWhite == subBoard->ScoreBestWhite && ioBoard->StepsBestWhite > subBoard->StepsBestWhite)) {
              ioBoard->StepsBestWhite= subBoard->StepsBestWhite + 1;
              ioBoard->ScoreBestWhite= subBoard->ScoreBestWhite;
              ioBoard->MoveBestWhite[0]= {w, h};
              ioBoard->MoveBestWhite[1]= ioBoard->Moves[w][h][k];
            }
          }
        }
        else {
          // TODO handle 2 player mode with minimax dependant on iDepth and idxTurn
        }
      }
    }
  }
}


void JumpinPlayerAI::ComputeBoardMoves(BoardState *ioBoard, const int iDepth) {
  if (ioBoard == nullptr) printf("Error: Computing moves of a null board");

  // Sweep through pawns
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      // Clear existing moves and boards
      ioBoard->Moves[w][h].clear();
      DeleteSubBoards(ioBoard, w, h);
      // Compute possible moves for the current pawn
      if (ioBoard->Pawns[w][h] != 0) {
        std::vector<std::vector<bool>> Visit= Field::AllocField2D(nW, nH, false);
        Visit[w][h]= true;
        const int oldPawn= ioBoard->Pawns[w][h];
        ioBoard->Pawns[w][h]= 0;
        ComputePawnDestinations(ioBoard, Visit, w, h, w, h);
        ioBoard->Pawns[w][h]= oldPawn;
      }
      // Create board for each move of the current pawn
      for (std::array<int, 2> move : ioBoard->Moves[w][h]) {
        BoardState *newBoard= CreateBoard();
        newBoard->Pawns= ioBoard->Pawns;
        newBoard->Pawns[move[0]][move[1]]= newBoard->Pawns[w][h];
        newBoard->Pawns[w][h]= 0;
        ioBoard->SubBoards[w][h].push_back(newBoard);
        if (iDepth < D.UI[SearchDepth_____].I() - 1) {
          ComputeBoardMoves(newBoard, iDepth + 1);
        }
      }
      // Compute score for each move of the current pawn
      for (BoardState *subBoard : ioBoard->SubBoards[w][h]) {
        ComputeBoardScore(subBoard);
      }
    }
  }
}


void JumpinPlayerAI::ComputePawnDestinations(BoardState *ioBoard,
                                             std::vector<std::vector<bool>> &ioVisit,
                                             const int iStartW, const int iStartH,
                                             const int iW, const int iH) {
  for (int idxDir= 0; idxDir < 4; idxDir++) {
    const int wInc= (idxDir == 0) ? (-1) : ((idxDir == 1) ? (+1) : (0));
    const int hInc= (idxDir == 2) ? (-1) : ((idxDir == 3) ? (+1) : (0));
    if (iW + 2 * wInc < 0 || iW + 2 * wInc >= nW ||
        iH + 2 * hInc < 0 || iH + 2 * hInc >= nH) continue;
    if (ioBoard->Pawns[iW + wInc][iH + hInc] == 0) continue;
    int wOff= iW + 2 * wInc;
    int hOff= iH + 2 * hInc;
    while (wOff >= 0 && wOff < nW && hOff >= 0 && hOff < nH) {
      if (ioBoard->Pawns[wOff][hOff] != 0) {
        wOff+= wInc;
        hOff+= hInc;
      }
      else if (!ioVisit[wOff][hOff]) {
        ioBoard->Moves[iStartW][iStartH].push_back(std::array<int, 2>{wOff, hOff});
        ioVisit[wOff][hOff]= true;
        ComputePawnDestinations(ioBoard, ioVisit, iStartW, iStartH, wOff, hOff);
        break;
      }
      else
        break;
    }
  }
}
