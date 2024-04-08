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
    D.UI.push_back(ParamUI("BoardW__________", 6));
    D.UI.push_back(ParamUI("BoardH__________", 11));
    D.UI.push_back(ParamUI("StartingRows____", 2));
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
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool JumpinPlayerAI::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void JumpinPlayerAI::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  nW= std::max(D.UI[BoardW__________].I(), 1);
  nH= std::max(D.UI[BoardH__________].I(), 1);

  const double voxSize= 1.0 / (double)nH;
  D.boxMin= {0.5 - 0.5 * voxSize,
             0.5 - 0.5 * voxSize * double(nW),
             0.5 - 0.5 * voxSize * double(nH)};
  D.boxMax= {0.5 + 0.5 * voxSize,
             0.5 + 0.5 * voxSize * double(nW),
             0.5 + 0.5 * voxSize * double(nH)};
}


// Refresh the project
void JumpinPlayerAI::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  Board.Pawns= Field::AllocField2D(nW, nH, char(0));
  Board.Moves= Field::AllocField2D(nW, nH, std::vector<std::array<char, 2>>());
  Board.Score= 0;

  Select= {0, 0};
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (h < D.UI[StartingRows____].I()) Board.Pawns[w][h]= -1;
      if (h >= nH - D.UI[StartingRows____].I()) Board.Pawns[w][h]= +1;
    }
  }
  ComputeScore();
  ComputeMoves();
}


// Handle keypress
void JumpinPlayerAI::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Get cursor coordinates on the board
  const int wCursor= std::min(std::max(int(std::floor((double)nW * (D.mouseProjX[1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]))), 0), nW - 1);
  const int hCursor= std::min(std::max(int(std::floor((double)nH * (D.mouseProjX[2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]))), 0), nH - 1);

  // Handle keyboard action
  if (key == 'E') {
    Board.Pawns[wCursor][hCursor]= 0;
    ComputeScore();
    ComputeMoves();
  }
  if (key == 'B') {
    Board.Pawns[wCursor][hCursor]= -1;
    ComputeScore();
    ComputeMoves();
  }
  if (key == 'W') {
    Board.Pawns[wCursor][hCursor]= +1;
    ComputeScore();
    ComputeMoves();
  }
  if (key == 'S') {
    Select= {char(wCursor), char(hCursor)};
  }
  if (key == 'D') {
    std::array<char, 2> target= {char(wCursor), char(hCursor)};
    for (std::array<char, 2> move : Board.Moves[Select[0]][Select[1]]) {
      if (move == target) {
        Board.Pawns[target[0]][target[1]]= Board.Pawns[Select[0]][Select[1]];
        Board.Pawns[Select[0]][Select[1]]= 0;
        ComputeScore();
        ComputeMoves();
        break;
      }
    }
  }
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

  const double voxSize= 1.0 / (double)nH;

  // Draw the board
  if (D.displayMode1) {
    glPushMatrix();
    glTranslatef(0.5f, 0.0f, 0.0f);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if ((w + h) % 2 == 0) glColor3f(0.4f, 0.4f, 0.4f);
        else glColor3f(0.6f, 0.6f, 0.6f);
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
        if (Board.Pawns[w][h] != 0) {
          if (Board.Pawns[w][h] < 0) glColor3f(0.2f, 0.2f, 0.2f);
          if (Board.Pawns[w][h] > 0) glColor3f(0.8f, 0.8f, 0.8f);
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

  // Draw the selection and available moves
  if (D.displayMode3) {
    glLineWidth(2.0);
    glPushMatrix();
    // Set the initial transformation
    glTranslatef(D.boxMin[0] + 0.5f * voxSize, D.boxMin[1] + 0.5f * voxSize, D.boxMin[2] + 0.5f * voxSize);
    glScalef(voxSize, voxSize, voxSize);
    // Draw the selection
    glColor3f(0.2f, 0.2f, 1.0f);
    glPushMatrix();
    glTranslatef(0.0f, (float)Select[0], (float)Select[1]);
    glutWireCube(0.95f);
    glPopMatrix();
    // Draw the moves
    glColor3f(1.0f, 0.2f, 0.2f);
    for (std::array<char, 2> move : Board.Moves[Select[0]][Select[1]]) {
      glPushMatrix();
      glTranslatef(0.0f, (float)move[0], (float)move[1]);
      glutWireCube(0.95f);
      glPopMatrix();
    }
    glPopMatrix();
    glLineWidth(1.0);
  }
}


void JumpinPlayerAI::ComputeScore() {
  Board.Score= 0;
}


void JumpinPlayerAI::ComputeMoves() {
  Board.Moves= Field::AllocField2D(nW, nH, std::vector<std::array<char, 2>>());
}


void JumpinPlayerAI::ComputeDestinations(const int w, const int h) {
  // Destinations[w][h]= true;
  // for (int idxDir= 0; idxDir < 4; idxDir++) {
  //   const int wInc= (idxDir == 0) ? (-1) : ((idxDir == 1) ? (+1) : (0));
  //   const int hInc= (idxDir == 2) ? (-1) : ((idxDir == 3) ? (+1) : (0));
  //   if (w + 2 * wInc < 0 || w + 2 * wInc >= nW) continue;
  //   if (h + 2 * hInc < 0 || h + 2 * hInc >= nH) continue;
  //   if (!Occupied[w + wInc][h + hInc]) continue;
  //   int wOff= w + 2 * wInc;
  //   int hOff= h + 2 * hInc;
  //   while (wOff >= 0 && wOff < nW && hOff >= 0 && hOff < nH) {
  //     if (!Occupied[wOff][hOff]) {
  //       if (!Destinations[wOff][hOff]) ComputeDestinations(wOff, hOff);
  //       break;
  //     }
  //     wOff+= wInc;
  //     hOff+= hInc;
  //   }
  // }
}
