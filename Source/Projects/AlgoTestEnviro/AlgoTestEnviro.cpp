#include "AlgoTestEnviro.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileOutput.hpp"
#include "Geom/BoxGrid.hpp"
#include "Geom/MarchingCubes.hpp"
#include "Geom/MergeVertices.hpp"
#include "Geom/PrimitiveCSG.hpp"
#include "Geom/Sketch.hpp"
#include "Math/Field.hpp"
#include "Math/Functions.hpp"
#include "Math/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
AlgoTestEnviro::AlgoTestEnviro() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void AlgoTestEnviro::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("TestParamALG_00_", 0.9));
    D.UI.push_back(ParamUI("TestParamALG_01_", 0.1));
    D.UI.push_back(ParamUI("TestParamALG_02_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_03_", 64.0));
    D.UI.push_back(ParamUI("TestParamALG_04_", 64.0));
    D.UI.push_back(ParamUI("TestParamALG_05_", 64.0));
    D.UI.push_back(ParamUI("TestParamALG_06_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_07_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_08_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_09_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_10_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_11_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_12_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_13_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_14_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_15_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_16_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_17_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_18_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_19_", 0.0));
    D.UI.push_back(ParamUI("Isocut__________", 0.0));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));
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
bool AlgoTestEnviro::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool AlgoTestEnviro::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void AlgoTestEnviro::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Allocate()\n");
}


// Refresh the project
void AlgoTestEnviro::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Refresh()\n");
}


// Handle keypress
void AlgoTestEnviro::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable

  if (D.UI[VerboseLevel____].I() >= 5) printf("KeyPress()\n");

  // Testing basis computation
  if (key == 'T') {
    Vec::Vec3<float> U(D.UI[TestParamALG_00_].D(), D.UI[TestParamALG_01_].D(), D.UI[TestParamALG_02_].D());
    Vec::Vec3<float> V, W;
    U.normalize();
    U.computeBasis(V, W);
    printf("U= %f %f %f\n", U[0], U[1], U[2]);
    printf("V= %f %f %f\n", V[0], V[1], V[2]);
    printf("Z= %f %f %f\n", W[0], W[1], W[2]);

    Verts.clear();
    Bars.clear();
    Tris.clear();

    Verts.push_back(std::array<double, 3>{0.0, 0.0, 0.0});
    Verts.push_back(std::array<double, 3>{U[0], U[1], U[2]});
    Verts.push_back(std::array<double, 3>{V[0], V[1], V[2]});
    Verts.push_back(std::array<double, 3>{W[0], W[1], W[2]});

    Bars.push_back(std::array<int, 2>{0, 1});
    Bars.push_back(std::array<int, 2>{0, 2});
    Bars.push_back(std::array<int, 2>{0, 3});

    Tris.push_back(std::array<int, 3>{0, 1, 2});
    Tris.push_back(std::array<int, 3>{0, 2, 3});
    Tris.push_back(std::array<int, 3>{0, 3, 1});
  }

  // Testing CSG field
  if (key == 'F') {
    ScalarField= Field::AllocField3D(D.UI[TestParamALG_03_].I(), D.UI[TestParamALG_04_].I(), D.UI[TestParamALG_05_].I(), std::numeric_limits<double>::max());
    std::array<double, 3> P0({0.1, 0.2, 0.3});
    std::array<double, 3> P1({0.8, 0.6, 0.7});
    PrimitiveCSG::ConeRound(P0, P1, 0.1, 0.2, PrimitiveCSG::BooleanMode::Union, D.boxMin, D.boxMax, ScalarField);
  }

  // Testing sketch smoothing
  if (key == 'S') {
    std::vector<std::array<double, 3>> polylineRef;
    polylineRef.push_back(std::array<double, 3>({0.0, 0.0, 0.0}));
    polylineRef.push_back(std::array<double, 3>({0.5, 0.5, 0.5}));
    polylineRef.push_back(std::array<double, 3>({1.0, -0.5, -0.5}));
    std::vector<std::array<double, 3>> polylineNew= polylineRef;
    Sketch::PolylineSubdivideAndSmooth(false, 3, 2, polylineNew);

    Verts.clear();
    Bars.clear();
    Tris.clear();
    for (int k= 0; k < (int)polylineNew.size(); k++) {
      Verts.push_back(polylineNew[k]);
      Bars.push_back(std::array<int, 2>{k, (k + 1) % (int)polylineNew.size()});
    }
  }

  // Testing Marching Cubes and vertex merge
  if (key == 'M') {
    Verts.clear();
    Bars.clear();
    Tris.clear();
    MarchingCubes::ComputeMarchingCubes(0.0, D.boxMin, D.boxMax, ScalarField, Verts, Tris);
    if (D.UI[TestParamALG_06_].B())
      MergeVertices::QuadraticMerge(D.UI[TestParamALG_07_].D(), Verts, Tris);
  }

  // Testing OBJ file output
  if (key == 'O') {
    std::vector<std::array<double, 3>> VertsCol;
    std::vector<std::array<int, 4>> Quads;
    FileOutput::SaveMeshOBJFile("FileOutput/test.obj", Verts, VertsCol, Tris, Quads, D.UI[VerboseLevel____].I());
  }

  // Testing Penal function
  if (key == 'P') {
    Verts.clear();
    for (int k= 0; k < 1000; k++) {
      const double ratio= double(k) / double(1000 - 1);
      Verts.push_back(std::array<double, 3>{D.boxMin[0] + (D.boxMax[0] - D.boxMin[0]) * 0.5,
                                            D.boxMin[1] + (D.boxMax[1] - D.boxMin[1]) * ratio,
                                            D.boxMin[2] + (D.boxMax[2] - D.boxMin[2]) * Functions::PenalInterpo(ratio, D.UI[TestParamALG_08_].D(), false)});
      Verts.push_back(std::array<double, 3>{D.boxMin[0] + (D.boxMax[0] - D.boxMin[0]) * 0.5,
                                            D.boxMin[1] + (D.boxMax[1] - D.boxMin[1]) * ratio,
                                            D.boxMin[2] + (D.boxMax[2] - D.boxMin[2]) * Functions::PenalInterpo(ratio, D.UI[TestParamALG_08_].D(), true)});
    }
  }
}


// Handle mouse action
void AlgoTestEnviro::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
}


// Animate the project
void AlgoTestEnviro::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
  if (D.UI[VerboseLevel____].I() >= 5) printf("Animate()\n");
}


// Draw the project
void AlgoTestEnviro::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;
  if (D.UI[VerboseLevel____].I() >= 5) printf("Draw()\n");

  // Draw vertices
  if (D.displayMode1) {
    glPointSize(5.0f);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_POINTS);
    for (int k= 0; k < (int)Verts.size(); k++) {
      glVertex3f(Verts[k][0], Verts[k][1], Verts[k][2]);
    }
    glEnd();
    glPointSize(1.0f);
  }

  // Draw bars
  if (D.displayMode2) {
    glColor3f(0.6f, 0.6f, 0.6f);
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    for (int k= 0; k < (int)Bars.size(); k++) {
      glVertex3f(Verts[Bars[k][0]][0], Verts[Bars[k][0]][1], Verts[Bars[k][0]][2]);
      glVertex3f(Verts[Bars[k][1]][0], Verts[Bars[k][1]][1], Verts[Bars[k][1]][2]);
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw triangles
  if (D.displayMode3) {
    glEnable(GL_LIGHTING);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_TRIANGLES);
    for (int k= 0; k < (int)Tris.size(); k++) {
      Vec::Vec3 v0(Verts[Tris[k][0]][0], Verts[Tris[k][0]][1], Verts[Tris[k][0]][2]);
      Vec::Vec3 v1(Verts[Tris[k][1]][0], Verts[Tris[k][1]][1], Verts[Tris[k][1]][2]);
      Vec::Vec3 v2(Verts[Tris[k][2]][0], Verts[Tris[k][2]][1], Verts[Tris[k][2]][2]);
      Vec::Vec3 n0= (v1 - v0).cross(v2 - v0).normalized();
      Vec::Vec3 n1= n0;
      Vec::Vec3 n2= n0;
      glNormal3f(n0[0], n0[1], n0[2]);
      glVertex3f(v0[0], v0[1], v0[2]);
      glNormal3f(n1[0], n1[1], n1[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
      glNormal3f(n2[0], n2[1], n2[2]);
      glVertex3f(v2[0], v2[1], v2[2]);
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }

  // Draw scalar field
  if (D.displayMode4) {
    int nbX, nbY, nbZ;
    Field::GetFieldDimensions(ScalarField, nbX, nbY, nbZ);
    if (nbX > 0 && nbY > 0 && nbZ > 0) {
      double stepX, stepY, stepZ, voxDiag, startX, startY, startZ;
      BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ, voxDiag);
      BoxGrid::GetVoxelStart(D.boxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      for (int x= 0; x < nbX; x++) {
        for (int y= 0; y < nbY; y++) {
          for (int z= 0; z < nbZ; z++) {
            if (ScalarField[x][y][z] < D.UI[Isocut__________].D()) continue;
            float r= 0.0f, g= 0.0f, b= 0.0f;
            Colormap::RatioToJetBrightSmooth(0.5 * (ScalarField[x][y][z] + 0.5) * D.UI[ColorFactor_____].D(), r, g, b);
            glColor3f(r, g, b);
            glVertex3f(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
          }
        }
      }
      glEnd();
      glPointSize(1.0f);
    }
  }
}
