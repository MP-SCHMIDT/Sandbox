#include "JumpinPlayerAI.hpp"


// Standard lib
#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
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
#include "Util/Timer.hpp"

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
    D.UI.push_back(ParamUI("BoardW__________", 6));
    D.UI.push_back(ParamUI("BoardH__________", 11));
    D.UI.push_back(ParamUI("StartingRows____", 2));
    D.UI.push_back(ParamUI("RandomBoard_____", 0));
    D.UI.push_back(ParamUI("BotBluPlayer____", 1));
    D.UI.push_back(ParamUI("BotRedPlayer____", 1));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("BluMoveOrdering_", 0));
    D.UI.push_back(ParamUI("RedMoveOrdering_", 1));
    D.UI.push_back(ParamUI("MaxSearchDepth__", 6));
    D.UI.push_back(ParamUI("MaxThinkTime____", 0.0));
    D.UI.push_back(ParamUI("MaxTreeBoards___", 8000));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("ValPushTotal____", 10));
    D.UI.push_back(ParamUI("ValPushLast_____", 50));
    D.UI.push_back(ParamUI("ValSoftStranded_", -40));
    D.UI.push_back(ParamUI("ValHardStranded_", -100));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.e-2));
    D.UI.push_back(ParamUI("______________03", NAN));
    D.UI.push_back(ParamUI("TestParamGAI_0__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_1__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_2__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_3__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_4__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_5__", 0.0));
    D.UI.push_back(ParamUI("______________04", NAN));
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
  if (D.UI[RandomBoard_____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool JumpinPlayerAI::CheckRefresh() {
  if (D.UI[BluMoveOrdering_].hasChanged()) isRefreshed= false;
  if (D.UI[RedMoveOrdering_].hasChanged()) isRefreshed= false;
  if (D.UI[MaxSearchDepth__].hasChanged()) isRefreshed= false;
  if (D.UI[MaxThinkTime____].hasChanged()) isRefreshed= false;
  if (D.UI[MaxTreeBoards___].hasChanged()) isRefreshed= false;
  if (D.UI[ValPushTotal____].hasChanged()) isRefreshed= false;
  if (D.UI[ValPushLast_____].hasChanged()) isRefreshed= false;
  if (D.UI[ValSoftStranded_].hasChanged()) isRefreshed= false;
  if (D.UI[ValHardStranded_].hasChanged()) isRefreshed= false;
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
  nbTreeBoards= 0;
  thinkTime= 0.0;
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
  std::vector<std::vector<int>> StartPawns= Field::AllocField2D(nW, nH, 0);
  int nbRows= std::min(D.UI[StartingRows____].I(), nH / 2);
  if (D.UI[RandomBoard_____].B()) {
    for (int player= -1; player <= +1; player+= 2) {
      for (int k= 0; k < nbRows * nW; k++) {
        while (true) {
          const int randW= Random::Val(0, nW - 1);
          const int randH= Random::Val(0, nH - 1);
          if (StartPawns[randW][randH] == 0) {
            StartPawns[randW][randH]= player;
            break;
          }
        }
      }
    }
  }
  else {
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (h < nbRows)
          StartPawns[w][h]= +1;
        if (h >= nH - nbRows)
          StartPawns[w][h]= -1;
      }
    }
  }
  RootBoard= CreateBoard(StartPawns, std::array<int, 4>({0, 0, 0, 0}));
}


// Refresh the project
void JumpinPlayerAI::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Compute root board score and moves
  ComputeGameTreeSearch();

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
      ComputeGameTreeSearch();
    }
  }

  // Skip turn
  if (key == 'C') {
    idxTurn++;
    ComputeGameTreeSearch();
  }

  // Manual pawn selection
  if (key == 'S') {
    wSel= wCursor;
    hSel= hCursor;
  }

  // Manual pawn move
  if (key == 'D') {
    if (wSel >= 0 && wSel < nW && hSel >= 0 && hSel < nH) {
      std::array<int, 4> manualMove= {wSel, hSel, wCursor, hCursor};
      for (BoardState *subBoard : RootBoard->SubBoards) {
        if (subBoard->Move == manualMove) {
          RootBoard->Pawns[wCursor][hCursor]= RootBoard->Pawns[wSel][hSel];
          RootBoard->Pawns[wSel][hSel]= 0;
          wSel= hSel= -1;
          idxTurn++;
          ComputeGameTreeSearch();
          break;
        }
      }
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

  // Autoplay the Nash move if the current player is a bot
  if ((D.UI[BotBluPlayer____].B() && IsBluTurn(0)) || (D.UI[BotRedPlayer____].B() && !IsBluTurn(0))) {
    if (!RootBoard->SubBoards.empty()) {
      std::array<int, 4> NashMove= RootBoard->SubBoards[0]->Move;
      RootBoard->Pawns[NashMove[2]][NashMove[3]]= RootBoard->Pawns[NashMove[0]][NashMove[1]];
      RootBoard->Pawns[NashMove[0]][NashMove[1]]= 0;
      wSel= hSel= -1;
      idxTurn++;
      ComputeGameTreeSearch();
    }
  }

  // Plot the scores
  PlotData();
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

  // TODO add cheat mode overlay showing Nash moves

  // Draw the available moves for the selected pawn
  if (D.displayMode3) {
    glLineWidth(2.0f);
    glPushMatrix();
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    if (wSel >= 0 && wSel < nW && hSel >= 0 && hSel < nH) {
      for (BoardState *subBoard : RootBoard->SubBoards) {
        if (subBoard->Move[0] == wSel && subBoard->Move[1] == hSel) {
          glColor3f(0.5f, 0.8f, 0.5f);
          glPushMatrix();
          glTranslatef(0.0f, (float)subBoard->Move[2], (float)subBoard->Move[3]);
          glutWireCube(0.95f);
          glPopMatrix();
        }
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
    DrawBoardTree(RootBoard, 0, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), -1.0f, std::numbers::pi + 1.0f);
    glEnd();
    glLineWidth(1.0f);
    glPointSize(10.0f);
    glBegin(GL_POINTS);
    DrawBoardTree(RootBoard, 0, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), -1.0f, std::numbers::pi + 1.0f);
    glEnd();
    glPointSize(1.0f);
  }
}


void JumpinPlayerAI::DrawBoardTree(const BoardState *iBoard, const int iDepth,
                                   const float px, const float py, const float pz,
                                   const float radius, const float arcBeg, const float arcEnd) {
  if (iBoard == nullptr) printf("[ERROR] DrawBoardTree on a null board\n");

  // Precompute distances and arc radians
  const float arcStep= (arcEnd - arcBeg) / float(iBoard->SubBoards.size());
  const float distBeg= radius * float(iDepth) / D.UI[MaxSearchDepth__].F();
  const float distEnd= radius * float(iDepth + 1) / D.UI[MaxSearchDepth__].F();

  // Recursively draw the moves
  for (int idxMove= 0; idxMove < (int)iBoard->SubBoards.size(); idxMove++) {
    float r, g, b;
    Colormap::RatioToJetBrightSmooth(0.5f - D.UI[ColorFactor_____].F() * float(iBoard->SubBoards[idxMove]->Score), r, g, b);
    glColor3f(r, g, b);
    glVertex3f(px, py + distBeg * std::cos((arcBeg + arcEnd) / 2.0f), pz + distBeg * std::sin((arcBeg + arcEnd) / 2.0f));
    glVertex3f(px, py + distEnd * std::cos(arcBeg + 0.5f * arcStep + idxMove * arcStep), pz + distEnd * std::sin(arcBeg + 0.5f * arcStep + idxMove * arcStep));
    DrawBoardTree(iBoard->SubBoards[idxMove], iDepth + 1, px, py, pz, radius, arcBeg + 0.5f * arcStep + (idxMove - 0.5f) * arcStep, arcBeg + 0.5f * arcStep + (idxMove + 0.5f) * arcStep);
  }
}


void JumpinPlayerAI::PlotData() {
  // Plot score evolution
  if (D.Plot.size() < 2) D.Plot.resize(2);
  D.Plot[0].name= "Score";
  D.Plot[0].isSymmetric= true;
  D.Plot[0].val.push_back(RootBoard->Score);
  D.Plot[1].name= "NashScore";
  D.Plot[1].isSymmetric= true;
  D.Plot[1].isSameRange= true;
  D.Plot[1].val.push_back(RootBoard->NashScore);

  // Print turn and win state
  D.Status.clear();
  D.Status.resize(4);
  D.Status[0]= std::string{"Turn:"} + std::to_string(idxTurn);
  D.Status[1]= std::string{"Player:"} + (IsBluTurn(0) ? std::string{"Blu"} : std::string{"Red"});
  if (RootBoard->NashScore == +INT_MAX) D.Status[2]= std::string{"BluWin:"} + std::to_string(RootBoard->NashNbSteps);
  if (RootBoard->NashScore == -INT_MAX) D.Status[2]= std::string{"RedWin:"} + std::to_string(RootBoard->NashNbSteps);
  D.Status[3]= std::string{"ThinkTime:"} + std::to_string(thinkTime) + std::string{"ms"};
}


JumpinPlayerAI::BoardState *JumpinPlayerAI::CreateBoard(const std::vector<std::vector<int>> &iPawns,
                                                        const std::array<int, 4> &iMove) {
  BoardState *newBoard= new BoardState;
  newBoard->Pawns= iPawns;
  newBoard->Move= iMove;
  newBoard->Score= 0;
  newBoard->NashScore= 0;
  newBoard->NashNbSteps= 0;
  return newBoard;
}


void JumpinPlayerAI::DeleteBoard(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] DeleteBoard on a null board\n");
  for (BoardState *subBoard : ioBoard->SubBoards)
    DeleteBoard(subBoard);
  delete ioBoard;
  ioBoard= nullptr;
}


void JumpinPlayerAI::ComputeBoardScore(BoardState *ioBoard) {
  if (ioBoard == nullptr) printf("[ERROR] ComputeBoardScore on a null board\n");

  // Check for win state
  bool isWinBlu= true;
  bool isWinRed= true;
  for (int w= 0; w < nW && isWinBlu; w++)
    for (int h= nH - D.UI[StartingRows____].I(); h < nH && isWinBlu; h++)
      if (ioBoard->Pawns[w][h] <= 0) isWinBlu= false;
  for (int w= 0; w < nW && isWinRed; w++)
    for (int h= 0; h < D.UI[StartingRows____].I() && isWinRed; h++)
      if (ioBoard->Pawns[w][h] >= 0) isWinRed= false;

  if (isWinBlu) ioBoard->Score= +INT_MAX;
  else if (isWinRed) ioBoard->Score= -INT_MAX;
  else {
    // Reset the score
    ioBoard->Score= 0;

    // Add score for all pawn advance
    if (D.UI[ValPushTotal____].I() != 0) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          if (ioBoard->Pawns[w][h] > 0) ioBoard->Score+= (h + 1) * D.UI[ValPushTotal____].I();
          if (ioBoard->Pawns[w][h] < 0) ioBoard->Score-= (nH - h) * D.UI[ValPushTotal____].I();
        }
      }
    }

    // Add score for last pawn advance
    if (D.UI[ValPushLast_____].I() != 0) {
      bool foundLastBlu= false;
      for (int h= 0; h < nH && !foundLastBlu; h++) {
        for (int w= 0; w < nW && !foundLastBlu; w++) {
          if (ioBoard->Pawns[w][h] > 0) {
            ioBoard->Score+= (h + 1) * D.UI[ValPushLast_____].I();
            foundLastBlu= true;
          }
        }
      }
      bool foundLastRed= false;
      for (int h= nH - 1; h >= 0 && !foundLastRed; h--) {
        for (int w= 0; w < nW && !foundLastRed; w++) {
          if (ioBoard->Pawns[w][h] < 0) {
            ioBoard->Score-= (nH - h) * D.UI[ValPushLast_____].I();
            foundLastRed= true;
          }
        }
      }
    }

    // Add score for soft stranded pawns
    if (D.UI[ValSoftStranded_].I() != 0) {
      // TODO
    }

    // Add score for hard stranded pawns
    if (D.UI[ValHardStranded_].I() != 0) {
      // TODO
    }
  }
}


void JumpinPlayerAI::ComputeGameTreeSearch() {
  // Start the timer and board counter
  Timer::PushTimer();
  nbTreeBoards= 1;

  // Reset root board
  for (BoardState *subBoard : RootBoard->SubBoards)
    DeleteBoard(subBoard);
  RootBoard->SubBoards.clear();
  RootBoard->Move= {0, 0, 0, 0};
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
  }
  // Handle first visit of board and creation of sub boards
  else if (ioBoard->SubBoards.empty()) {
    // Find all the possible moves from the current board
    std::vector<std::array<int, 4>> Moves;
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if ((ioBoard->Pawns[w][h] > 0 && IsBluTurn(iDepth)) || (ioBoard->Pawns[w][h] < 0 && !IsBluTurn(iDepth))) {
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
    // Create and recursively search the sub board for each move
    for (std::array<int, 4> move : Moves) {
      // Create and score the sub board for the current move
      BoardState *newBoard= CreateBoard(ioBoard->Pawns, move);
      newBoard->Pawns[move[2]][move[3]]= newBoard->Pawns[move[0]][move[1]];
      newBoard->Pawns[move[0]][move[1]]= 0;
      ComputeBoardScore(newBoard);
      ioBoard->SubBoards.push_back(newBoard);
      nbTreeBoards++;
      // Recursively search the sub tree for the current move
      RecursiveTreeSearch(newBoard, iDepth + 1, iMaxDepth);
    }
    // Sort the sub boards
    if (!ioBoard->SubBoards.empty()) {
      SortSubBoards(ioBoard, iDepth);
      ioBoard->NashScore= ioBoard->SubBoards[0]->NashScore;
      ioBoard->NashNbSteps= ioBoard->SubBoards[0]->NashNbSteps + 1;
    }
  }
  // Handle already visited board with existing sub boards
  else {
    for (BoardState *subBoard : ioBoard->SubBoards)
      RecursiveTreeSearch(subBoard, iDepth + 1, iMaxDepth);
    SortSubBoards(ioBoard, iDepth);
    ioBoard->NashScore= ioBoard->SubBoards[0]->NashScore;
    ioBoard->NashNbSteps= ioBoard->SubBoards[0]->NashNbSteps + 1;
  }
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

  // Get the ordering mode
  int ordering= 0;
  if (IsBluTurn(iDepth)) ordering= D.UI[BluMoveOrdering_].I();
  else ordering= D.UI[RedMoveOrdering_].I();

  // Bubble sort of the sub boards in Nash order
  if (ordering == 0) {
    for (int k0= 0; k0 < (int)ioBoard->SubBoards.size(); k0++) {
      for (int k1= k0 + 1; k1 < (int)ioBoard->SubBoards.size(); k1++) {
        if ((IsBluTurn(iDepth) && ioBoard->SubBoards[k0]->NashScore < ioBoard->SubBoards[k1]->NashScore) ||
            (!IsBluTurn(iDepth) && ioBoard->SubBoards[k0]->NashScore > ioBoard->SubBoards[k1]->NashScore) ||
            (ioBoard->SubBoards[k0]->NashScore == ioBoard->SubBoards[k1]->NashScore &&
             ioBoard->SubBoards[k0]->NashNbSteps > ioBoard->SubBoards[k1]->NashNbSteps)) {
          BoardState *tmp= ioBoard->SubBoards[k0];
          ioBoard->SubBoards[k0]= ioBoard->SubBoards[k1];
          ioBoard->SubBoards[k1]= tmp;
        }
      }
    }
  }

  // Bubble sort of the sub boards in greedy order
  else if (ordering == 1) {
    for (int k0= 0; k0 < (int)ioBoard->SubBoards.size(); k0++) {
      for (int k1= k0 + 1; k1 < (int)ioBoard->SubBoards.size(); k1++) {
        if ((IsBluTurn(iDepth) && ioBoard->SubBoards[k0]->Score < ioBoard->SubBoards[k1]->Score) ||
            (!IsBluTurn(iDepth) && ioBoard->SubBoards[k0]->Score > ioBoard->SubBoards[k1]->Score)) {
          BoardState *tmp= ioBoard->SubBoards[k0];
          ioBoard->SubBoards[k0]= ioBoard->SubBoards[k1];
          ioBoard->SubBoards[k1]= tmp;
        }
      }
    }
  }

  // Random shuffle of the boards
  else if (ordering == 2) {
    static auto rng= std::default_random_engine{};
    std::shuffle(std::begin(ioBoard->SubBoards), std::end(ioBoard->SubBoards), rng);
  }
}


bool inline JumpinPlayerAI::IsBluTurn(const int iDepth) {
  return (idxTurn + iDepth) % 2 == 0;
}
