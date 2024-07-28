#include "TheBoardGameAI.hpp"


// Standard lib
#include <cmath>
#include <format>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
TheBoardGameAI::TheBoardGameAI() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void TheBoardGameAI::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("GameMode________", 0));
    D.UI.push_back(ParamUI("BoardW__________", 6));
    D.UI.push_back(ParamUI("BoardH__________", 6));
    D.UI.push_back(ParamUI("MoveStreakRed___", 1));
    D.UI.push_back(ParamUI("MoveStreakBlu___", 1));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("MaxSearchDepth__", 4));
    D.UI.push_back(ParamUI("MaxThinkTime____", 0.0));
    D.UI.push_back(ParamUI("MaxTreeBoards___", 0));
    D.UI.push_back(ParamUI("MoveSortScore___", 1));
    D.UI.push_back(ParamUI("MoveSortNash____", 1));
    D.UI.push_back(ParamUI("MoveSortRand____", 1));
    D.UI.push_back(ParamUI("ABPruning_______", 1));
    D.UI.push_back(ParamUI("IterDeepening___", 1));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("HexEdgeConnect__", 10));
    D.UI.push_back(ParamUI("HexCornConnect__", 25));
    D.UI.push_back(ParamUI("JmpPushTotal____", 1));
    D.UI.push_back(ParamUI("JmpPushLast_____", 0));
    D.UI.push_back(ParamUI("JmpSoftStranded_", 0));
    D.UI.push_back(ParamUI("JmpHardStranded_", 0));
    D.UI.push_back(ParamUI("ChkMaterial_____", 10));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("RandomMoves_____", 0));
    D.UI.push_back(ParamUI("BotStrategyRed__", 3));
    D.UI.push_back(ParamUI("BotStrategyBlu__", 3));
    D.UI.push_back(ParamUI("______________03", NAN));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.e-2));
    D.UI.push_back(ParamUI("______________04", NAN));
    D.UI.push_back(ParamUI("TestParamGAI_0__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_1__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_2__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_3__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_4__", 0.0));
    D.UI.push_back(ParamUI("TestParamGAI_5__", 0.0));
    D.UI.push_back(ParamUI("______________05", NAN));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));

    D.displayModeLabel[1]= "Board";
    D.displayModeLabel[2]= "Pawns";
    D.displayModeLabel[3]= "Moves";
    D.displayModeLabel[4]= "Tree";
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool TheBoardGameAI::CheckAlloc() {
  if (D.UI[GameMode________].hasChanged()) isAllocated= false;
  if (D.UI[BoardW__________].hasChanged()) isAllocated= false;
  if (D.UI[BoardH__________].hasChanged()) isAllocated= false;
  if (D.UI[MoveStreakRed___].hasChanged()) isAllocated= false;
  if (D.UI[MoveStreakBlu___].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool TheBoardGameAI::CheckRefresh() {
  if (D.UI[MaxSearchDepth__].hasChanged()) isRefreshed= false;
  if (D.UI[MaxThinkTime____].hasChanged()) isRefreshed= false;
  if (D.UI[MaxTreeBoards___].hasChanged()) isRefreshed= false;
  if (D.UI[MoveSortScore___].hasChanged()) isRefreshed= false;
  if (D.UI[MoveSortNash____].hasChanged()) isRefreshed= false;
  if (D.UI[MoveSortRand____].hasChanged()) isRefreshed= false;
  if (D.UI[ABPruning_______].hasChanged()) isRefreshed= false;
  if (D.UI[IterDeepening___].hasChanged()) isRefreshed= false;

  if (D.UI[HexEdgeConnect__].hasChanged()) isRefreshed= false;
  if (D.UI[HexCornConnect__].hasChanged()) isRefreshed= false;
  if (D.UI[JmpPushTotal____].hasChanged()) isRefreshed= false;
  if (D.UI[JmpPushLast_____].hasChanged()) isRefreshed= false;
  if (D.UI[JmpSoftStranded_].hasChanged()) isRefreshed= false;
  if (D.UI[JmpHardStranded_].hasChanged()) isRefreshed= false;
  if (D.UI[ChkMaterial_____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void TheBoardGameAI::Allocate() {
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
  streakRed= std::max(D.UI[MoveStreakRed___].I(), 0);
  streakBlu= std::max(D.UI[MoveStreakBlu___].I(), 0);
  if (streakBlu == 0) streakRed= std::max(streakRed, 1);
  nbTreeBoards= 0;
  thinkTime= 0.0;
  wSel= -1;
  hSel= -1;

  D.boxMin= {0.5, 0.0, 0.0};
  D.boxMax= {0.5, 1.0, 1.0};

  // Build the board cell positions
  if (D.UI[GameMode________].I() == 0) {
    const float wStep= 1.0f;
    const float hStep= std::sin(std::numbers::pi / 3.0f);
    const float cloudWidth= float(nW - 1) * wStep + float(nH - 1) * 0.5f * hStep + wStep;
    const float cloudHeight= float(nH - 1) * hStep + hStep;
    cellSize= 1.0f / std::max(cloudWidth, cloudHeight);
    Cells= Field::Field2(nW, nH, Vec::Vec3<float>(0.0, 0.0, 0.0));
    for (int w= 0; w < nW; w++)
      for (int h= 0; h < nH; h++)
        Cells.at(w, h).set(0.5f,
                           0.5f + ((float(w) + 0.5f) * wStep + float(nH - 1 - h) * 0.5f * wStep - 0.5f * cloudWidth) * cellSize,
                           0.5f - ((float(nH - 1 - h) + 0.5f) * hStep - 0.5f * cloudHeight) * cellSize);
  }
  else if (D.UI[GameMode________].I() == 1) {
    cellSize= 1.0f / (float)std::max(nW, nH);
    Cells= Field::Field2(nW, nH, Vec::Vec3<float>(0.0, 0.0, 0.0));
    for (int w= 0; w < nW; w++)
      for (int h= 0; h < nH; h++)
        Cells.at(w, h).set(0.5f,
                           0.5f * cellSize + float(w) * cellSize,
                           0.5f * cellSize + float(h) * cellSize);
  }
  else if (D.UI[GameMode________].I() == 2) {
    cellSize= 1.0f / (float)std::max(nW, nH);
    Cells= Field::Field2(nW, nH, Vec::Vec3<float>(0.0, 0.0, 0.0));
    for (int w= 0; w < nW; w++)
      for (int h= 0; h < nH; h++)
        Cells.at(w, h).set(0.5f,
                           0.5f * cellSize + float(w) * cellSize,
                           0.5f * cellSize + float(h) * cellSize);
  }

  // Initialize the pawns
  Field::Field2<int> Pawns(nW, nH, 0);
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (D.UI[GameMode________].I() == 0) {
      }
      else if (D.UI[GameMode________].I() == 1) {
        if (streakRed > 0 && h < 2) Pawns.at(w, h)= +1;
        if (streakBlu > 0 && h >= nH - 2) Pawns.at(w, h)= -1;
      }
      else if (D.UI[GameMode________].I() == 2) {
        if ((w + h) % 2 == 0) {
          if (h < 3) Pawns.at(w, h)= +1;
          if (h >= nH - 3) Pawns.at(w, h)= -1;
        }
      }
    }
  }

  // Create the root board
  RootBoard= CreateBoard(Pawns, std::vector<std::array<int, 2>>(), 0);
  if (D.UI[GameMode________].I() == 0) ComputeBoardScoreHex(RootBoard);
  if (D.UI[GameMode________].I() == 1) ComputeBoardScoreJmp(RootBoard);
  if (D.UI[GameMode________].I() == 2) ComputeBoardScoreChk(RootBoard);
  RootBoard->NashScore= RootBoard->Score;
  RootBoard->NashNbSteps= 0;
}


// Refresh the project
void TheBoardGameAI::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Compute root board score and moves
  ComputeGameTreeSearch(D.UI[MaxSearchDepth__].I());

  // Plot the scores
  PlotData();
}


// Handle UI parameter change
void TheBoardGameAI::ParamChange() {
}


// Handle keypress
void TheBoardGameAI::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Get cursor coordinates on the board
  int wCursor= 0;
  int hCursor= 0;
  Vec::Vec3<float> mouseProj(D.mouseProjX[0], D.mouseProjX[1], D.mouseProjX[2]);
  float distMin= (mouseProj - Cells.at(0, 0)).norm();
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      const float dist= (mouseProj - Cells.at(w, h)).norm();
      if (distMin > dist) {
        distMin= dist;
        wCursor= w;
        hCursor= h;
      }
    }
  }

  // Manual board set-up
  if (D.keyLetterUpperCase == 'R' || D.keyLetterUpperCase == 'G' || D.keyLetterUpperCase == 'B') {
    if (D.keyLetterUpperCase == 'R') RootBoard->Pawns.at(wCursor, hCursor)= +1;
    if (D.keyLetterUpperCase == 'G') RootBoard->Pawns.at(wCursor, hCursor)= 0;
    if (D.keyLetterUpperCase == 'B') RootBoard->Pawns.at(wCursor, hCursor)= -1;
    wSel= hSel= -1;
    ComputeGameTreeSearch(D.UI[MaxSearchDepth__].I());
  }

  // Skip turn
  if (D.keyLetterUpperCase == 'C') {
    idxTurn++;
    wSel= hSel= -1;
    ComputeGameTreeSearch(D.UI[MaxSearchDepth__].I());
  }

  // Play a sequence random moves
  if (D.keyLetterUpperCase == 'A') {
    for (int k= 0; k < D.UI[RandomMoves_____].I(); k++) {
      ComputeGameTreeSearch(1);
      if (!RootBoard->SubBoards.empty()) {
        const int idxMove= Random::Val(0, (int)RootBoard->SubBoards.size() - 1);
        RootBoard->Move= RootBoard->SubBoards[idxMove]->Move;
        if (D.UI[GameMode________].I() == 0) ExecuteMoveHex(RootBoard, 0);
        else if (D.UI[GameMode________].I() == 1) ExecuteMoveJmp(RootBoard, 0);
        else if (D.UI[GameMode________].I() == 2) ExecuteMoveChk(RootBoard, 0);
      }
      idxTurn++;
    }
    ComputeGameTreeSearch(D.UI[MaxSearchDepth__].I());
  }

  // Run benchmark of current AI capability
  if (D.keyLetterUpperCase == 'T') {
    // Start with win position
    // Execute N random moves
    // Run AI to see if it returns to win position in N or more moves
  }

  // Manual pawn selection
  if (D.keyLetterUpperCase == 'S') {
    wSel= wCursor;
    hSel= hCursor;
  }

  // Manual pawn placement
  if (D.keyLetterUpperCase == 'D') {
    // Find and execute the first move that matches the selection
    for (BoardState *subBoard : RootBoard->SubBoards) {
      if (D.UI[GameMode________].I() == 0) {
        wSel= subBoard->Move[0][0];
        hSel= subBoard->Move[0][1];
      }
      if (subBoard->Move[0] == std::array<int, 2>{wSel, hSel} &&
          subBoard->Move[subBoard->Move.size() - 1] == std::array<int, 2>{wCursor, hCursor}) {
        RootBoard->Move= subBoard->Move;
        if (D.UI[GameMode________].I() == 0) ExecuteMoveHex(RootBoard, 0);
        else if (D.UI[GameMode________].I() == 1) ExecuteMoveJmp(RootBoard, 0);
        else if (D.UI[GameMode________].I() == 2) ExecuteMoveChk(RootBoard, 0);
        wSel= hSel= -1;
        idxTurn++;
        ComputeGameTreeSearch(D.UI[MaxSearchDepth__].I());
        break;
      }
    }
  }

  // Plot the scores
  PlotData();
}


// Handle mouse action
void TheBoardGameAI::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void TheBoardGameAI::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Autoplay a move if the current player is a bot
  const int BotStrategy= IsRedTurn(0) ? D.UI[BotStrategyRed__].I() : D.UI[BotStrategyBlu__].I();
  if (BotStrategy > 0) {
    if (!RootBoard->SubBoards.empty()) {
      // Select the move based on chosen strategy
      // TODO add strategy for mostly nash move with fraction of random moves
      int idxMove= 0;
      if (BotStrategy == 1) {
        idxMove= Random::Val(0, (int)RootBoard->SubBoards.size() - 1);
      }
      else if (BotStrategy == 2) {
        idxMove= GetIdxBestSubBoard(RootBoard, 0, 0);
      }
      else if (BotStrategy == 3) {
        idxMove= GetIdxBestSubBoard(RootBoard, 0, 1);
      }
      // Execute the move
      RootBoard->Move= RootBoard->SubBoards[idxMove]->Move;
      if (D.UI[GameMode________].I() == 0) {
        ExecuteMoveHex(RootBoard, 0);
        ComputeBoardScoreHex(RootBoard);
      }
      if (D.UI[GameMode________].I() == 1) {
        ExecuteMoveJmp(RootBoard, 0);
        ComputeBoardScoreJmp(RootBoard);
      }
      if (D.UI[GameMode________].I() == 2) {
        ExecuteMoveChk(RootBoard, 0);
        ComputeBoardScoreChk(RootBoard);
      }
    }
    idxTurn++;
    ComputeGameTreeSearch(D.UI[MaxSearchDepth__].I());
  }

  // Plot the scores
  PlotData();
}


// Draw the project
void TheBoardGameAI::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the board
  if (D.displayMode[1]) {
    if (D.UI[GameMode________].I() == 0) {
      glLineWidth(8.0f);
      glBegin(GL_LINES);
      glColor3f(1.0f, 0.0f, 0.0f);
      glVertex3fv((Cells.at(0, 0) - Vec::Vec3<float>(0.0f, 0.0f, 0.7f * cellSize)).array());
      glVertex3fv((Cells.at(nW - 1, 0) - Vec::Vec3<float>(0.0f, 0.0f, 0.7f * cellSize)).array());
      glVertex3fv((Cells.at(0, nH - 1) + Vec::Vec3<float>(0.0f, 0.0f, 0.7f * cellSize)).array());
      glVertex3fv((Cells.at(nW - 1, nH - 1) + Vec::Vec3<float>(0.0f, 0.0f, 0.7f * cellSize)).array());
      glColor3f(0.0f, 0.0f, 1.0f);
      glVertex3fv((Cells.at(0, 0) - Vec::Vec3<float>(0.0f, 0.7f * cellSize, 0.0f)).array());
      glVertex3fv((Cells.at(0, nH - 1) - Vec::Vec3<float>(0.0f, 0.7f * cellSize, 0.0f)).array());
      glVertex3fv((Cells.at(nW - 1, 0) + Vec::Vec3<float>(0.0f, 0.7f * cellSize, 0.0f)).array());
      glVertex3fv((Cells.at(nW - 1, nH - 1) + Vec::Vec3<float>(0.0f, 0.7f * cellSize, 0.0f)).array());
      glEnd();
      glLineWidth(1.0f);
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          glColor3f(0.5f, 0.5f, 0.5f);
          glPushMatrix();
          glTranslatef(Cells.at(w, h)[0], Cells.at(w, h)[1], Cells.at(w, h)[2]);
          glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
          glScalef(1.0f, 1.0f, 0.01f);
          glutSolidSphere(0.5 * cellSize, 6, 2);
          glPopMatrix();
        }
      }
    }
    else if (D.UI[GameMode________].I() == 1) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          const float col= ((w + h) % 2 == 0) ? 0.4f : 0.6f;
          if (h < 2) glColor3f(col, col, col + 0.2f);
          else if (h >= nH - 2) glColor3f(col + 0.2f, col, col);
          else glColor3f(col, col, col);
          glPushMatrix();
          glTranslatef(Cells.at(w, h)[0], Cells.at(w, h)[1], Cells.at(w, h)[2]);
          glScalef(0.01f, 1.0f, 1.0f);
          glutSolidCube(cellSize);
          glPopMatrix();
        }
      }
    }
    else if (D.UI[GameMode________].I() == 2) {
      for (int w= 0; w < nW; w++) {
        for (int h= 0; h < nH; h++) {
          const float col= ((w + h) % 2 == 0) ? 0.4f : 0.6f;
          glColor3f(col, col, col);
          glPushMatrix();
          glTranslatef(Cells.at(w, h)[0], Cells.at(w, h)[1], Cells.at(w, h)[2]);
          glScalef(0.01f, 1.0f, 1.0f);
          glutSolidCube(cellSize);
          glPopMatrix();
        }
      }
    }
  }

  // Draw the pawns
  if (D.displayMode[2]) {
    glEnable(GL_LIGHTING);
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (RootBoard->Pawns.at(w, h) == 0) continue;
        if (RootBoard->Pawns.at(w, h) > 0) glColor3f(1.0f, 0.5f, 0.5f);
        if (RootBoard->Pawns.at(w, h) < 0) glColor3f(0.5f, 0.5f, 1.0f);
        glPushMatrix();
        glTranslatef(Cells.at(w, h)[0], Cells.at(w, h)[1], Cells.at(w, h)[2]);
        glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
        glutSolidSphere(0.4 * cellSize, 36, 10);
        glPopMatrix();
      }
    }
    glDisable(GL_LIGHTING);
  }

  // Draw the possible moves
  if (D.displayMode[3]) {
    for (BoardState *subBoard : RootBoard->SubBoards) {
      float r, g, b;
      if (D.UI[ColorMode_______].I() == 0) Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(subBoard->Score), r, g, b);
      if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(subBoard->NashScore), r, g, b);
      if (D.UI[ColorMode_______].I() == 2) Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(subBoard->NashNbSteps), r, g, b);
      glColor3f(r, g, b);
      const int w= subBoard->Move[subBoard->Move.size() - 1][0];
      const int h= subBoard->Move[subBoard->Move.size() - 1][1];
      glPushMatrix();
      glTranslatef(Cells.at(w, h)[0], Cells.at(w, h)[1], Cells.at(w, h)[2]);
      glutWireCube(0.5 * cellSize);
      glPopMatrix();
    }
  }

  // Draw the board tree
  if (D.displayMode[4]) {
    float px= 0.5f * (D.boxMin[0] + D.boxMax[0]);
    float py= D.boxMax[1] + 0.1f * (D.boxMax[1] - D.boxMin[1]);
    float pz= D.boxMin[2] + 0.5f * (D.boxMax[2] - D.boxMin[2]);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    DrawBoardTree(RootBoard, 0, 0, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), -0.5f * std::numbers::pi, 0.5f * std::numbers::pi);
    glEnd();
    glLineWidth(1.0f);
    glPointSize(6.0f);
    glBegin(GL_POINTS);
    DrawBoardTree(RootBoard, 0, 1, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), -0.5f * std::numbers::pi, 0.5f * std::numbers::pi);
    glEnd();
    glPointSize(1.0f);
  }
}


void TheBoardGameAI::DrawBoardTree(const BoardState *iBoard, const int iDepth, const int iDrawMode,
                                   const float px, const float py, const float pz,
                                   const float radius, const float arcBeg, const float arcEnd) {
  if (iBoard == nullptr) printf("[ERROR] DrawBoardTree on a null board\n");

  // Precompute distances and arc radians
  const float arcStep= (arcEnd - arcBeg) / float(iBoard->SubBoards.size());
  const float distBeg= radius * float(iDepth) / D.UI[MaxSearchDepth__].F();
  const float distEnd= radius * float(iDepth + 1) / D.UI[MaxSearchDepth__].F();

  // Recursively draw the moves
  for (int idxMove= 0; idxMove < (int)iBoard->SubBoards.size(); idxMove++) {
    if (iDrawMode == 0) {
      float r, g, b;
      if (D.UI[ColorMode_______].I() == 0) Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(iBoard->SubBoards[idxMove]->Score), r, g, b);
      if (D.UI[ColorMode_______].I() == 1) Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(iBoard->SubBoards[idxMove]->NashScore), r, g, b);
      if (D.UI[ColorMode_______].I() == 2) Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(iBoard->SubBoards[idxMove]->NashNbSteps), r, g, b);
      glColor3f(r, g, b);
      glVertex3f(px, py + distBeg * std::cos((arcBeg + arcEnd) / 2.0f), pz + distBeg * std::sin((arcBeg + arcEnd) / 2.0f));
      glVertex3f(px, py + distEnd * std::cos(arcBeg + 0.5f * arcStep + idxMove * arcStep), pz + distEnd * std::sin(arcBeg + 0.5f * arcStep + idxMove * arcStep));
    }
    else {
      if (IsRedTurn(iDepth)) glColor3f(1.0f, 0.0f, 0.0f);
      else glColor3f(0.0f, 0.0f, 1.0f);
      glVertex3f(px, py + distBeg * std::cos((arcBeg + arcEnd) / 2.0f), pz + distBeg * std::sin((arcBeg + arcEnd) / 2.0f));
    }
    DrawBoardTree(iBoard->SubBoards[idxMove], iDepth + 1, iDrawMode, px, py, pz, radius, arcBeg + 0.5f * arcStep + (idxMove - 0.5f) * arcStep, arcBeg + 0.5f * arcStep + (idxMove + 0.5f) * arcStep);
  }
}


void TheBoardGameAI::PlotData() {
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
  D.Status.resize(5);
  D.Status[0]= std::format("Turn:{}", idxTurn);
  D.Status[1]= std::format("ThinkTime:{:6.3f}ms", thinkTime);
  D.Status[2]= std::format("BoardCount:{}", nbTreeBoards);
  D.Status[3]= std::format("Player:{}", (IsRedTurn(0) ? 'R' : 'B'));
  if (RootBoard->NashScore == +INT_MAX) D.Status[4]= std::format("RedWin:{}", RootBoard->NashNbSteps);
  if (RootBoard->NashScore == -INT_MAX) D.Status[4]= std::format("BluWin:{}", RootBoard->NashNbSteps);
}
