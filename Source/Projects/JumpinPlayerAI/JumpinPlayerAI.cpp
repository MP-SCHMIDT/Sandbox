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
    D.UI.push_back(ParamUI("EnableBluPlayer_", 1));
    D.UI.push_back(ParamUI("EnableRedPlayer_", 1));
    D.UI.push_back(ParamUI("StartingRows____", 2));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("SearchDepth_____", 3));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.e-2));
    D.UI.push_back(ParamUI("______________02", NAN));
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
  if (D.UI[EnableBluPlayer_].hasChanged()) isAllocated= false;
  if (D.UI[EnableRedPlayer_].hasChanged()) isAllocated= false;
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
      if (D.UI[EnableBluPlayer_].B() && h < D.UI[StartingRows____].I())
        RootBoard->Pawns[w][h]= +1;
      if (D.UI[EnableRedPlayer_].B() && h >= nH - D.UI[StartingRows____].I())
        RootBoard->Pawns[w][h]= -1;
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
  ComputeGameTreeSearch(RootBoard, 0);

  // Plot the scores
  PlotData();
}


// Handle keypress
void JumpinPlayerAI::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Get cursor coordinates on the board
  const int wCursor= (int)std::floor((double)nW * (D.mouseProjX[1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
  const int hCursor= (int)std::floor((double)nH * (D.mouseProjX[2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));

  // Board set-up
  if (key == 'R' || key == 'G' || key == 'B') {
    if (wCursor >= 0 && wCursor < nW && hCursor >= 0 && hCursor < nH) {
      if (key == 'R') RootBoard->Pawns[wCursor][hCursor]= -1;
      if (key == 'G') RootBoard->Pawns[wCursor][hCursor]= 0;
      if (key == 'B') RootBoard->Pawns[wCursor][hCursor]= +1;
      wSel= hSel= -1;
      ComputeGameTreeSearch(RootBoard, 0);
    }
  }

  // Pawn selection
  if (key == 'S') {
    wSel= wCursor;
    hSel= hCursor;
  }

  // Skip turn
  if (key == 'C') {
    idxTurn++;
    ComputeGameTreeSearch(RootBoard, 0);
  }

  // Manual move
  bool AutoplayBot= false;
  if (key == 'D') {
    std::array<int, 2> target= {wCursor, hCursor};
    if (wSel >= 0 && wSel < nW && hSel >= 0 && hSel < nH) {
      for (std::array<int, 2> move : RootBoard->Moves[wSel][hSel]) {
        if (move == target) {
          RootBoard->Pawns[target[0]][target[1]]= RootBoard->Pawns[wSel][hSel];
          RootBoard->Pawns[wSel][hSel]= 0;
          wSel= hSel= -1;
          idxTurn++;
          ComputeGameTreeSearch(RootBoard, 0);
          AutoplayBot= true;
          break;
        }
      }
    }
  }

  // Automated move
  if (key == 'F' || AutoplayBot) {
    AutoplayBot= false;
    if (RootBoard->BestIsSet && RootBoard->StepsBest > 0) {
      std::array<int, 2> beg= {-1, -1}, end= {-1, -1};
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          for (int k= 0; k < (int)RootBoard->Moves[w][h].size(); k++) {
            BoardState *subBoard= RootBoard->SubBoards[w][h][k];
            if (RootBoard->ScoreBest == subBoard->ScoreBest && RootBoard->StepsBest == subBoard->StepsBest + 1) {
              beg= {w, h};
              end= RootBoard->Moves[w][h][k];
            }
          }
        }
      }
      RootBoard->Pawns[end[0]][end[1]]= RootBoard->Pawns[beg[0]][beg[1]];
      RootBoard->Pawns[beg[0]][beg[1]]= 0;
      wSel= hSel= -1;
      idxTurn++;
      ComputeGameTreeSearch(RootBoard, 0);
    }
  }

  // Plot the scores
  PlotData();
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
        const float col= ((w + h) % 2 == 0) ? 0.4f : 0.6f;
        if (h < D.UI[StartingRows____].I()) glColor3f(col + 0.2f, col, col);
        else if (h >= nH - D.UI[StartingRows____].I()) glColor3f(col, col, col + 0.2f);
        else glColor3f(col, col, col);
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
          float r= 0.5f, g= 0.5f, b= 0.5f;
          if (RootBoard->Pawns[w][h] < 0) r+= 0.4f;
          if (RootBoard->Pawns[w][h] > 0) b+= 0.4f;
          glColor3f(r, g, b);
          glPushMatrix();
          glTranslatef(D.boxMin[0] + 0.5f * voxSize,
                       D.boxMin[1] + 0.5f * voxSize + float(w) * voxSize,
                       D.boxMin[2] + 0.5f * voxSize + float(h) * voxSize);
          glutSolidSphere(0.45 * voxSize, 36, 10);
          glPopMatrix();
          if (w == wSel && h == hSel) {
            glColor3f(1.0f, 1.0f, 1.0f);
            glPushMatrix();
            glTranslatef(D.boxMin[0] + 0.9f * voxSize,
                         D.boxMin[1] + 0.5f * voxSize + float(w) * voxSize,
                         D.boxMin[2] + 0.5f * voxSize + float(h) * voxSize);
            glutSolidSphere(0.20 * voxSize, 36, 10);
            glPopMatrix();
          }
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
      for (std::array<int, 2> move : RootBoard->Moves[wSel][hSel]) {
        glColor3f(0.5f, 0.8f, 0.5f);
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
    DrawBoardTree(RootBoard, 0, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), 0.0f, 2.0f * std::numbers::pi);
    glEnd();
    glLineWidth(1.0f);
  }
}


void JumpinPlayerAI::DrawBoardTree(const BoardState *iBoard, const int iDepth,
                                   const float px, const float py, const float pz,
                                   const float radius, const float arcBeg, const float arcEnd) {
  if (iBoard == nullptr) printf("Error: Drawing a null board");

  // Find the number of moves for angle spacing
  int nbMoves= 0;
  for (int w= 0; w < nW; w++)
    for (int h= 0; h < nH; h++)
      nbMoves+= (int)iBoard->Moves[w][h].size();
  const float arcStep= (arcEnd - arcBeg) / float(nbMoves);
  const float distBeg= radius * float(iDepth) / D.UI[SearchDepth_____].F();
  const float distEnd= radius * float(iDepth + 1) / D.UI[SearchDepth_____].F();

  // Recursively draw the moves
  int idxMove= 0;
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      for (int k= 0; k < (int)iBoard->Moves[w][h].size(); k++) {
        BoardState *subBoard= iBoard->SubBoards[w][h][k];
        // Draw the branch for the move
        float r, g, b;
        Colormap::RatioToJetBrightSmooth(0.5f - D.UI[ColorFactor_____].F() * float(subBoard->Score), r, g, b);
        glColor3f(r, g, b);
        glVertex3f(px, py + distBeg * std::cos((arcBeg + arcEnd) / 2.0f), pz + distBeg * std::sin((arcBeg + arcEnd) / 2.0f));
        glVertex3f(px, py + distEnd * std::cos(arcBeg + idxMove * arcStep), pz + distEnd * std::sin(arcBeg + idxMove * arcStep));
        // Recursively draw moves in the sub arc
        DrawBoardTree(subBoard, iDepth + 1, px, py, pz, radius, arcBeg + (idxMove - 0.5f) * arcStep, arcBeg + (idxMove + 0.5f) * arcStep);
        idxMove++;
      }
    }
  }
}


void JumpinPlayerAI::PlotData() {
  if (D.Plot.size() < 2) D.Plot.resize(2);
  D.Plot[0].name= "Score";
  D.Plot[0].isSymmetric= true;
  D.Plot[0].val.push_back(RootBoard->Score);
  D.Plot[1].name= "BestScore";
  D.Plot[1].isSymmetric= true;
  D.Plot[1].isSameRange= true;
  D.Plot[1].val.push_back(RootBoard->ScoreBest);

  bool isBluTurn= true;
  if (D.UI[EnableRedPlayer_].B() && (!D.UI[EnableBluPlayer_].B() || idxTurn % 2 == 1)) isBluTurn= false;
  if ((int)D.Status.size() != 2) D.Status.resize(2);
  D.Status[0]= std::string{"PlayerTurn: "} + (isBluTurn ? std::string{"Blu"} : std::string{"Red"});
  D.Status[1]= std::string{"WinState: No"};
  if (RootBoard->WinState != 0) {
    D.Status[1]= std::string{"WinState: "} +
                 (RootBoard->WinState > 0 ? std::string{"Blu"} : std::string{"Red"}) +
                 std::string{" in "} + std::to_string(RootBoard->StepsBest);
  }
}


JumpinPlayerAI::BoardState *JumpinPlayerAI::CreateBoard() {
  BoardState *newBoard= new BoardState;
  newBoard->Pawns= Field::AllocField2D(nW, nH, 0);
  newBoard->Moves= Field::AllocField2D(nW, nH, std::vector<std::array<int, 2>>());
  newBoard->SubBoards= Field::AllocField2D(nW, nH, std::vector<BoardState *>());
  newBoard->Score= 0;
  newBoard->ScoreBest= 0;
  newBoard->StepsBest= 0;
  newBoard->WinState= 0;
  newBoard->BestIsSet= false;
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


void JumpinPlayerAI::ComputeGameTreeSearch(BoardState *ioBoard, const int iDepth) {
  if (ioBoard == nullptr) printf("Error: Computing game tree search of a null board");

  // Reset if root node
  if (iDepth == 0) {
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        ioBoard->Moves[w][h].clear();
        for (BoardState *subBoard : ioBoard->SubBoards[w][h])
          DeleteBoard(subBoard);
        ioBoard->SubBoards[w][h].clear();
      }
    }
    ioBoard->Score= 0;
    ioBoard->ScoreBest= 0;
    ioBoard->StepsBest= 0;
    ioBoard->WinState= 0;
    ioBoard->BestIsSet= false;
  }

  // Set the base score of the board
  ComputeBoardScore(ioBoard);
  CheckWinState(ioBoard, iDepth);
  if (iDepth == D.UI[SearchDepth_____].I() || ioBoard->WinState != 0) {
    ioBoard->ScoreBest= ioBoard->Score;
    ioBoard->StepsBest= 0;
    ioBoard->BestIsSet= true;
    return;
  }

  // Search the tree recursively
  ioBoard->BestIsSet= false;
  bool isBluTurn= true;
  if (D.UI[EnableRedPlayer_].B() && (!D.UI[EnableBluPlayer_].B() || (idxTurn + iDepth) % 2 == 1)) isBluTurn= false;
  if (iDepth < D.UI[SearchDepth_____].I()) {
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        // Compute possible moves for the current pawn
        if (ioBoard->Pawns[w][h] != 0) {
          if ((isBluTurn && ioBoard->Pawns[w][h] > 0) ||
              (!isBluTurn && ioBoard->Pawns[w][h] < 0)) {
            std::vector<std::vector<bool>> Visit= Field::AllocField2D(nW, nH, false);
            Visit[w][h]= true;
            const int oldPawn= ioBoard->Pawns[w][h];
            ioBoard->Pawns[w][h]= 0;
            ComputePawnJumps(ioBoard, Visit, w, h, w, h);
            ioBoard->Pawns[w][h]= oldPawn;
          }
        }
        // Create board for each move of the current pawn
        for (std::array<int, 2> move : ioBoard->Moves[w][h]) {
          BoardState *newBoard= CreateBoard();
          newBoard->Pawns= ioBoard->Pawns;
          newBoard->Pawns[move[0]][move[1]]= newBoard->Pawns[w][h];
          newBoard->Pawns[w][h]= 0;
          ioBoard->SubBoards[w][h].push_back(newBoard);
          if (iDepth < D.UI[SearchDepth_____].I()) {
            ComputeGameTreeSearch(newBoard, iDepth + 1);
          }
          // Update the best score
          if (isBluTurn) {
            if (!ioBoard->BestIsSet ||
                ioBoard->ScoreBest < newBoard->ScoreBest ||
                (ioBoard->ScoreBest == newBoard->ScoreBest && ioBoard->StepsBest > newBoard->StepsBest + 1)) {
              ioBoard->ScoreBest= newBoard->ScoreBest;
              ioBoard->StepsBest= newBoard->StepsBest + 1;
              ioBoard->WinState= newBoard->WinState;
              ioBoard->BestIsSet= true;
            }
          }
          else {
            if (!ioBoard->BestIsSet ||
                ioBoard->ScoreBest > newBoard->ScoreBest ||
                (ioBoard->ScoreBest == newBoard->ScoreBest && ioBoard->StepsBest > newBoard->StepsBest + 1)) {
              ioBoard->ScoreBest= newBoard->ScoreBest;
              ioBoard->StepsBest= newBoard->StepsBest + 1;
              ioBoard->WinState= newBoard->WinState;
              ioBoard->BestIsSet= true;
            }
          }
        }
      }
    }
  }
}


void JumpinPlayerAI::CheckWinState(BoardState *ioBoard, const int iDepth) {
  if (ioBoard == nullptr) printf("Error: Checking win state of a null board");

  bool isBluTurn= true;
  if (D.UI[EnableRedPlayer_].B() && (!D.UI[EnableBluPlayer_].B() || (idxTurn + iDepth) % 2 == 1)) isBluTurn= false;

  if (isBluTurn) {
    bool win= true;
    for (int w= 0; w < nW && win; w++)
      for (int h= nH - D.UI[StartingRows____].I(); h < nH && win; h++)
        if (ioBoard->Pawns[w][h] <= 0) win= false;
    if (win)
      ioBoard->WinState= +1;
  }
  else {
    bool win= true;
    for (int w= 0; w < nW && win; w++)
      for (int h= 0; h < D.UI[StartingRows____].I() && win; h++)
        if (ioBoard->Pawns[w][h] >= 0) win= false;
    if (win)
      ioBoard->WinState= -1;
  }
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

  // TODO add other scoring schemes
}


void JumpinPlayerAI::ComputePawnJumps(BoardState *ioBoard,
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
      if (ioBoard->Pawns[wOff][hOff] == 0) {
        if (!ioVisit[wOff][hOff]) {
          ioBoard->Moves[iStartW][iStartH].push_back(std::array<int, 2>{wOff, hOff});
          ioVisit[wOff][hOff]= true;
          ComputePawnJumps(ioBoard, ioVisit, iStartW, iStartH, wOff, hOff);
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
