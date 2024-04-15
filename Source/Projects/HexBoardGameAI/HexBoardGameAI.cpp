#include "HexBoardGameAI.hpp"


// Standard lib
#include <cmath>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
HexBoardGameAI::HexBoardGameAI() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void HexBoardGameAI::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("BoardW__________", 5));
    D.UI.push_back(ParamUI("BoardH__________", 5));
    D.UI.push_back(ParamUI("RandomBoard_____", 0));
    D.UI.push_back(ParamUI("SinglePlayer____", 0));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("BotStrategyRed__", 3));
    D.UI.push_back(ParamUI("BotStrategyBlu__", 3));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("MaxSearchDepth__", 4));
    D.UI.push_back(ParamUI("MaxThinkTime____", 0.0));
    D.UI.push_back(ParamUI("MaxTreeBoards___", 0));
    D.UI.push_back(ParamUI("MoveSortScore___", 0));
    D.UI.push_back(ParamUI("MoveSortNash____", 0));
    D.UI.push_back(ParamUI("MoveSortRand____", 0));
    D.UI.push_back(ParamUI("ABPruning_______", 0));
    D.UI.push_back(ParamUI("IterDeepening___", 0));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("ValEdgeConnect__", 1));
    D.UI.push_back(ParamUI("ValCornConnect__", 2));
    D.UI.push_back(ParamUI("______________03", NAN));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.e-2));
    D.UI.push_back(ParamUI("______________04", NAN));
    D.UI.push_back(ParamUI("TestParamHEX_0__", 0.0));
    D.UI.push_back(ParamUI("TestParamHEX_1__", 0.0));
    D.UI.push_back(ParamUI("TestParamHEX_2__", 0.0));
    D.UI.push_back(ParamUI("TestParamHEX_3__", 0.0));
    D.UI.push_back(ParamUI("TestParamHEX_4__", 0.0));
    D.UI.push_back(ParamUI("TestParamHEX_5__", 0.0));
    D.UI.push_back(ParamUI("______________05", NAN));
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
bool HexBoardGameAI::CheckAlloc() {
  if (D.UI[BoardW__________].hasChanged()) isAllocated= false;
  if (D.UI[BoardH__________].hasChanged()) isAllocated= false;
  if (D.UI[RandomBoard_____].hasChanged()) isAllocated= false;
  if (D.UI[SinglePlayer____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool HexBoardGameAI::CheckRefresh() {
  if (D.UI[MaxSearchDepth__].hasChanged()) isRefreshed= false;
  if (D.UI[MaxThinkTime____].hasChanged()) isRefreshed= false;
  if (D.UI[MaxTreeBoards___].hasChanged()) isRefreshed= false;
  if (D.UI[MoveSortScore___].hasChanged()) isRefreshed= false;
  if (D.UI[MoveSortNash____].hasChanged()) isRefreshed= false;
  if (D.UI[MoveSortRand____].hasChanged()) isRefreshed= false;
  if (D.UI[ABPruning_______].hasChanged()) isRefreshed= false;
  if (D.UI[IterDeepening___].hasChanged()) isRefreshed= false;
  if (D.UI[ValEdgeConnect__].hasChanged()) isRefreshed= false;
  if (D.UI[ValCornConnect__].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void HexBoardGameAI::Allocate() {
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

  // Build the Hex board positions
  Cells= Field::AllocField2D(nW, nH, Vec::Vec3<float>(0.0, 0.0, 0.0));
  for (int w= 0; w < nW; w++)
    for (int h= 0; h < nH; h++)
      Cells[w][h].set(0.0f, float(w) * sqrt(3.0f) - float(h) * 0.5f * sqrt(3.0f), float(h) * 3.0f / 2.0f);
  Vec::Vec3<float> boxMin= Cells[0][0];
  Vec::Vec3<float> boxMax= Cells[0][0];
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      boxMin= boxMin.cwiseMin(Cells[w][h]);
      boxMax= boxMax.cwiseMax(Cells[w][h]);
    }
  }
  cellSize= 1.0f / (boxMax - boxMin).maxCoeff();
  for (int w= 0; w < nW; w++)
    for (int h= 0; h < nH; h++)
      Cells[w][h]= Vec::Vec3<float>(0.f, 0.5f, 0.5f) + ((Cells[w][h] - 0.5f * (boxMax + boxMin)) * cellSize);
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      boxMin= boxMin.cwiseMin(Cells[w][h]);
      boxMax= boxMax.cwiseMax(Cells[w][h]);
    }
  }

  D.boxMin= {0.0, 0.0, 0.0};
  D.boxMax= {1.0, 1.0, 1.0};

  // Create and initialize root board with pawns
  std::vector<std::vector<int>> Pawns= Field::AllocField2D(nW, nH, 0);
  if (D.UI[RandomBoard_____].B()) {
    for (int player= -1; player <= +1; player+= 2) {
      for (int k= 0; k < D.UI[RandomBoard_____].I(); k++) {
        while (true) {
          const int randW= Random::Val(0, nW - 1);
          const int randH= Random::Val(0, nH - 1);
          if (Pawns[randW][randH] == 0) {
            Pawns[randW][randH]= player;
            break;
          }
        }
      }
    }
  }

  // Remove Blu for single player
  if (D.UI[SinglePlayer____].B()) {
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (Pawns[w][h] == -1)
          Pawns[w][h]= 0;
      }
    }
  }

  RootBoard= CreateBoard(Pawns, std::array<int, 2>({-1, -1}), 0);
  ComputeBoardScore(RootBoard);
  RootBoard->NashScore= RootBoard->Score;
  RootBoard->NashNbSteps= 0;
}


// Refresh the project
void HexBoardGameAI::Refresh() {
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
void HexBoardGameAI::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Get cursor coordinates on the board
  int wCursor= 0;
  int hCursor= 0;
  Vec::Vec3<float> mouseProj(D.mouseProjX[0], D.mouseProjX[1], D.mouseProjX[2]);
  float distMin= (mouseProj - Cells[0][0]).norm();
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      const float dist= (mouseProj - Cells[w][h]).norm();
      if (distMin > dist) {
        distMin= dist;
        wCursor= w;
        hCursor= h;
      }
    }
  }

  // Board set-up
  if (key == 'R' || key == 'G' || key == 'B') {
    if (key == 'R') RootBoard->Pawns[wCursor][hCursor]= +1;
    if (key == 'G') RootBoard->Pawns[wCursor][hCursor]= 0;
    if (key == 'B') RootBoard->Pawns[wCursor][hCursor]= -1;
    ComputeGameTreeSearch();
  }

  // Skip turn
  if (key == 'C') {
    idxTurn++;
    ComputeGameTreeSearch();
  }

  // Run benchmark of current AI capability
  if (key == 'T') {
    // Start with win position
    // Execute N random moves
    // Run AI to see if it returns to win position in N or more moves
  }

  // Manual pawn placement
  if (key == 'D') {
    if (RootBoard->Pawns[wCursor][hCursor] == 0) {
      if (IsRedTurn(0)) RootBoard->Pawns[wCursor][hCursor]= +1;
      else RootBoard->Pawns[wCursor][hCursor]= -1;
      idxTurn++;
      ComputeGameTreeSearch();
    }
  }

  // Plot the scores
  PlotData();
}


// Animate the project
void HexBoardGameAI::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Autoplay a move if the current player is a bot
  const int BotStrategy= IsRedTurn(0) ? D.UI[BotStrategyRed__].I() : D.UI[BotStrategyBlu__].I();
  if (BotStrategy > 0) {
    if (!RootBoard->SubBoards.empty()) {
      // Select the move based on chosen strategy
      // TODO add stratyegy for mostly nash move with fraction of random moves
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
      std::array<int, 2> SelectedMove= RootBoard->SubBoards[idxMove]->Move;
      if (SelectedMove[0] != -1 && SelectedMove[1] != -1) {
        if (IsRedTurn(0)) RootBoard->Pawns[SelectedMove[0]][SelectedMove[1]]= +1;
        else RootBoard->Pawns[SelectedMove[0]][SelectedMove[1]]= -1;
      }
    }
    idxTurn++;
    ComputeGameTreeSearch();
  }

  // Plot the scores
  PlotData();
}


// Draw the project
void HexBoardGameAI::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the pawns
  if (D.displayMode1) {
    glEnable(GL_LIGHTING);
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        float r= 0.5f, g= 0.5f, b= 0.5f;
        if (RootBoard->Pawns[w][h] > 0) r+= 0.4f;
        if (RootBoard->Pawns[w][h] < 0) b+= 0.4f;
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(Cells[w][h][0], Cells[w][h][1], Cells[w][h][2]);
        glutSolidSphere(0.45 * cellSize, 36, 10);
        glPopMatrix();
      }
    }
    glDisable(GL_LIGHTING);
  }

  // TODO add cheat mode overlay showing Nash moves

  // Draw the board tree
  if (D.displayMode4) {
    float px= 0.5f * (D.boxMin[0] + D.boxMax[0]);
    float py= D.boxMax[1] + 0.6f * (D.boxMax[1] - D.boxMin[1]);
    float pz= D.boxMin[2] + 0.5f * (D.boxMax[2] - D.boxMin[2]);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    DrawBoardTree(RootBoard, 0, 0, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), -std::numbers::pi, 0.0f);
    glEnd();
    glLineWidth(1.0f);
    glPointSize(6.0f);
    glBegin(GL_POINTS);
    DrawBoardTree(RootBoard, 0, 1, px, py, pz, 0.5f * (D.boxMax[1] - D.boxMin[1]), -std::numbers::pi, 0.0f);
    glEnd();
    glPointSize(1.0f);
  }
}


void HexBoardGameAI::DrawBoardTree(const BoardState *iBoard, const int iDepth, const int iDrawMode,
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
      Colormap::RatioToJetBrightSmooth(0.5f + D.UI[ColorFactor_____].F() * float(iBoard->SubBoards[idxMove]->Score), r, g, b);
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


void HexBoardGameAI::PlotData() {
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
  D.Status[0]= std::string{"Turn:"} + std::to_string(idxTurn);
  D.Status[1]= std::string{"ThinkTime:"} + std::to_string(thinkTime) + std::string{"ms"};
  D.Status[2]= std::string{"BoardCount:"} + std::to_string(nbTreeBoards);
  D.Status[3]= std::string{"Player:"} + (IsRedTurn(0) ? std::string{"Red"} : std::string{"Blu"});
  if (RootBoard->NashScore == +INT_MAX) D.Status[4]= std::string{"RedWin:"} + std::to_string(RootBoard->NashNbSteps);
  if (RootBoard->NashScore == -INT_MAX) D.Status[4]= std::string{"BluWin:"} + std::to_string(RootBoard->NashNbSteps);
}
