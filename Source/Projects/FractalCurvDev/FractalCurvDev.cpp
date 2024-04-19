#include "FractalCurvDev.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Math/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
FractalCurvDev::FractalCurvDev() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void FractalCurvDev::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("FractalMode_____", 0));
    D.UI.push_back(ParamUI("MaxDepth________", 5));
    D.UI.push_back(ParamUI("StepZVal________", 0.2));
    D.UI.push_back(ParamUI("StepZExpo_______", 1.5));
    D.UI.push_back(ParamUI("SpreadCoeff_____", 1.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));
  }

  if (D.UI.size() != VerboseLevel____ + 1) {
    printf("[ERROR] Invalid parameter count in UI\n");
  }

  D.boxMin= {0.0, 0.0, 0.0};
  D.boxMax= {1.0, 1.0, 1.0};

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool FractalCurvDev::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool FractalCurvDev::CheckRefresh() {
  if (D.UI[FractalMode_____].hasChanged()) isRefreshed= false;
  if (D.UI[MaxDepth________].hasChanged()) isRefreshed= false;
  if (D.UI[StepZVal________].hasChanged()) isRefreshed= false;
  if (D.UI[StepZExpo_______].hasChanged()) isRefreshed= false;
  if (D.UI[SpreadCoeff_____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void FractalCurvDev::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;
}


// Refresh the project
void FractalCurvDev::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  int maxDepth= std::max(D.UI[MaxDepth________].I(), 2);
  Faces.clear();
  Nodes.clear();

  // Initialize the curve at depth 0
  if (D.UI[FractalMode_____].I() == 0) {
    // Koch snowflake
    Nodes.resize(1);
    Nodes[0].push_back(Vec::Vec3<float>(-std::sqrt(3.0f) / 2.0f, -1.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(+std::sqrt(3.0f) / 2.0f, +0.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(-std::sqrt(3.0f) / 2.0f, +1.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(-std::sqrt(3.0f) / 2.0f, -1.0f, 0.0f));
  }
  if (D.UI[FractalMode_____].I() == 1) {
    // Dragon curve
    Nodes.resize(1);
    Nodes[0].push_back(Vec::Vec3<float>(+0.0f, -1.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(+1.0f, +0.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(+0.0f, +1.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(-1.0f, +0.0f, 0.0f));
    Nodes[0].push_back(Vec::Vec3<float>(+0.0f, -1.0f, 0.0f));
  }

  // Iteratively build the next level in the fractal recursion depth
  for (int idxDepth= 1; idxDepth < maxDepth; idxDepth++) {
    Nodes.push_back(std::vector<Vec::Vec3<float>>());
    for (int idxNode= 0; idxNode < int(Nodes[idxDepth - 1].size()) - 1; idxNode++) {
      Vec::Vec3<float> posA= Nodes[idxDepth - 1][idxNode];
      Vec::Vec3<float> posB= Nodes[idxDepth - 1][idxNode + 1];
      Vec::Vec3<float> ZOffset(0.0f, 0.0f, -D.UI[StepZVal________].F() / std::pow(D.UI[StepZExpo_______].F(), double(idxDepth)));

      if (D.UI[FractalMode_____].I() == 0) {
        // Koch snowflake
        Vec::Vec3<float> posN0= ZOffset + posA;
        Vec::Vec3<float> posN1= ZOffset + posA + (posB - posA) * 1.0 / 3.0;
        Vec::Vec3<float> posN2= ZOffset + (posA + posB) / 2.0 + (std::sqrt(3.0) / 6.0) * (posB - posA).norm() * (posB - posA).cross(Vec::Vec3<float>(0.0, 0.0, 1.0)).normalized();
        Vec::Vec3<float> posN3= ZOffset + posA + (posB - posA) * 2.0 / 3.0;
        Vec::Vec3<float> posN4= ZOffset + posB;

        if (idxNode == 0)
          Nodes[idxDepth].push_back(posN0);
        Nodes[idxDepth].push_back(posN1);
        Nodes[idxDepth].push_back(posN2);
        Nodes[idxDepth].push_back(posN3);
        Nodes[idxDepth].push_back(posN4);

        Vec::Vec3<float> posM= (posA + posB) / 2.0f;
        Faces.push_back({posM, posA, posN0});
        Faces.push_back({posM, posN0, posN1});
        Faces.push_back({posM, posN1, posN2});
        Faces.push_back({posM, posN2, posN3});
        Faces.push_back({posM, posN3, posN4});
        Faces.push_back({posM, posN4, posB});
      }

      if (D.UI[FractalMode_____].I() == 1) {
        // Dragon curve
        Vec::Vec3<float> dir= (posB - posA).cross(Vec::Vec3<float>(0.0f, 0.0f, 1.0f)).normalized();
        if (idxNode % 2 == 0)
          dir= Vec::Vec3<float>() - dir;
        Vec::Vec3<float> posN0= ZOffset + posA;
        Vec::Vec3<float> posN1= ZOffset + 0.5f * (posA + posB) + 0.5f * D.UI[SpreadCoeff_____].F() * (posB - posA).norm() * dir;
        Vec::Vec3<float> posN2= ZOffset + posB;

        if (idxNode == 0)
          Nodes[idxDepth].push_back(posN0);
        Nodes[idxDepth].push_back(posN1);
        Nodes[idxDepth].push_back(posN2);

        Vec::Vec3<float> posM= (posA + posB) / 2.0f;
        Faces.push_back({posM, posA, posN0});
        Faces.push_back({posM, posN0, posN1});
        Faces.push_back({posM, posN1, posN2});
        Faces.push_back({posM, posN2, posB});
      }
    }
  }
}


// Handle keypress
void FractalCurvDev::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Save OBJ file of developed surface
  if (key == 'O') {
    std::string iFullpath= "FileOutput/test.obj";
    printf("Saving OBJ mesh file [%s]\n", iFullpath.c_str());
    FILE* outputFile= nullptr;
    outputFile= fopen(iFullpath.c_str(), "w");
    if (outputFile == nullptr) {
      printf("[ERROR] Unable to create the file\n\n");
    }
    else {
      for (unsigned int k= 0; k < Faces.size(); k++) {
        fprintf(outputFile, "v %lf %lf %lf\n", Faces[k][0][0], Faces[k][0][1], Faces[k][0][2]);
        fprintf(outputFile, "v %lf %lf %lf\n", Faces[k][1][0], Faces[k][1][1], Faces[k][1][2]);
        fprintf(outputFile, "v %lf %lf %lf\n", Faces[k][2][0], Faces[k][2][1], Faces[k][2][2]);
      }
      for (unsigned int k= 0; k < Faces.size(); k++) {
        fprintf(outputFile, "f %d %d %d\n", k * 3 + 1, k * 3 + 2, k * 3 + 3);
      }
      fclose(outputFile);
      printf("File saved: %zd vertices, %zd triangles\n", Faces.size() * 3, Faces.size());
    }
  }
}


// Handle mouse action
void FractalCurvDev::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
}


// Animate the project
void FractalCurvDev::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
}


// Draw the project
void FractalCurvDev::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw vertices
  if (D.displayMode1) {
    glPointSize(2.0f);
    glBegin(GL_POINTS);
    for (int idxDepth= 0; idxDepth < int(Nodes.size()); idxDepth++) {
      float r, g, b;
      Colormap::RatioToJetBrightSmooth(float(idxDepth) / float(Nodes.size() - 1), r, g, b);
      glColor3f(r, g, b);
      for (int idxNode= 0; idxNode < int(Nodes[idxDepth].size()); idxNode++) {
        glVertex3fv(Nodes[idxDepth][idxNode].array());
      }
    }
    glEnd();
    glPointSize(1.0f);
  }

  // Draw wireframe
  if (D.displayMode2) {
    glLineWidth(2.0f);
    for (int idxDepth= 0; idxDepth < int(Nodes.size()); idxDepth++) {
      float r, g, b;
      Colormap::RatioToJetBrightSmooth(float(idxDepth) / float(Nodes.size() - 1), r, g, b);
      glBegin(GL_LINE_STRIP);
      glColor3f(r, g, b);
      for (int idxNode= 0; idxNode < int(Nodes[idxDepth].size()); idxNode++) {
        glVertex3fv(Nodes[idxDepth][idxNode].array());
      }
      glEnd();
    }
    glLineWidth(1.0f);
  }

  // Draw faces
  if (D.displayMode3) {
    glEnable(GL_LIGHTING);
    glColor3f(0.3f, 0.3f, 0.3f);
    glBegin(GL_TRIANGLES);
    for (auto face : Faces) {
      Vec::Vec3<float> normal= (face[1] - face[0]).cross(face[2] - face[0]).normalized();
      glNormal3fv(normal.array());
      glVertex3fv(face[0].array());
      glVertex3fv(face[1].array());
      glVertex3fv(face[2].array());
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }
}
