#include "JumpinPlayerAI.hpp"


// Standard lib
#include <cmath>
#include <cstdlib>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Field.hpp"
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
    D.UI.push_back(ParamUI("BoardH__________", 8));
    D.UI.push_back(ParamUI("StartingRows____", 2));
    D.UI.push_back(ParamUI("RandomBoard_____", 0));
    D.UI.push_back(ParamUI("SinglePlayer____", 1));
    D.UI.push_back(ParamUI("BotStrategyBlu__", 3));
    D.UI.push_back(ParamUI("BotStrategyRed__", 3));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("MaxSearchDepth__", 6));
    D.UI.push_back(ParamUI("MaxThinkTime____", 0.0));
    D.UI.push_back(ParamUI("MaxTreeBoards___", 64000));
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
  if (D.UI[SinglePlayer____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool JumpinPlayerAI::CheckRefresh() {
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
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (h < nbRows) StartPawns[w][h]= +1;
      if (h >= nH - nbRows) StartPawns[w][h]= -1;
    }
  }

  // Overwrite with random positions
  if (D.UI[RandomBoard_____].B()) {
    srand(0);
    StartPawns= Field::AllocField2D(nW, nH, 0);
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

  // Remove red for single player
  if (D.UI[SinglePlayer____].B()) {
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (StartPawns[w][h] == -1)
          StartPawns[w][h]= 0;
      }
    }
  }

  RootBoard= CreateBoard(StartPawns, std::array<int, 4>({0, 0, 0, 0}), 0);
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

  // Run benchmark of current AI capability
  if (key == 'T') {
    // Start with win position
    // Execute N random moves
    // Run AI to see if it returns to win position in N or more moves
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

  // Autoplay a move if the current player is a bot
  const int BotStrategy= IsBluTurn(0) ? D.UI[BotStrategyBlu__].I() : D.UI[BotStrategyRed__].I();
  if (BotStrategy > 0) {
    if (!RootBoard->SubBoards.empty()) {
      // Select the move based on chosen strategy
      int idxMove= 0;
      if (BotStrategy == 1) {
        idxMove= Random::Val(0, (int)RootBoard->SubBoards.size() - 1);
      }
      else if (BotStrategy == 2) {
        for (int k= 0; k < (int)RootBoard->SubBoards.size(); k++) {
          if (IsBluTurn(0) && RootBoard->SubBoards[k]->Score > RootBoard->SubBoards[idxMove]->Score) idxMove= k;
          if (!IsBluTurn(0) && RootBoard->SubBoards[k]->Score < RootBoard->SubBoards[idxMove]->Score) idxMove= k;
        }
      }
      else if (BotStrategy == 3) {
        idxMove= 0;
      }
      // Execute the move
      std::array<int, 4> SelectedMove= RootBoard->SubBoards[idxMove]->Move;
      if (SelectedMove[0] != -1 && SelectedMove[1] != -1 && SelectedMove[2] != -1 && SelectedMove[3] != -1) {
        RootBoard->Pawns[SelectedMove[2]][SelectedMove[3]]= RootBoard->Pawns[SelectedMove[0]][SelectedMove[1]];
        RootBoard->Pawns[SelectedMove[0]][SelectedMove[1]]= 0;
      }
    }

    wSel= hSel= -1;
    idxTurn++;
    ComputeGameTreeSearch();
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
  D.Status[1]= std::string{"ThinkTime:"} + std::to_string(thinkTime) + std::string{"ms"};
  D.Status[2]= std::string{"Player:"} + (IsBluTurn(0) ? std::string{"Blu"} : std::string{"Red"});
  if (RootBoard->NashScore == +INT_MAX) D.Status[3]= std::string{"BluWin:"} + std::to_string(RootBoard->NashNbSteps);
  if (RootBoard->NashScore == -INT_MAX) D.Status[3]= std::string{"RedWin:"} + std::to_string(RootBoard->NashNbSteps);
}
