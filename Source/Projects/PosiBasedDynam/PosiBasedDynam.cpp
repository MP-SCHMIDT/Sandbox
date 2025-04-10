#include "PosiBasedDynam.hpp"


// Standard lib
#include <cmath>
#include <vector>
#include <random>
#include <map>

// Algo headers
#include "Type/Vec.hpp"
#include "FileIO/FileInput.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
PosiBasedDynam::PosiBasedDynam() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void PosiBasedDynam::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("DomainX_________", 0.5));
    D.UI.push_back(ParamUI("DomainY_________", 1.4));
    D.UI.push_back(ParamUI("DomainZ_________", 1.5));
    D.UI.push_back(ParamUI("NumParticl______", 1024));
    D.UI.push_back(ParamUI("RadParticl______", 0.03));
    D.UI.push_back(ParamUI("RadMouse________", 0.25));
    D.UI.push_back(ParamUI("AddTetMeshID____", 2));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("TimeStep________", 0.01));
    D.UI.push_back(ParamUI("NbSubSteps______", 4));
    D.UI.push_back(ParamUI("MaterialDensity_", 1.0));
    D.UI.push_back(ParamUI("ForceDrag_______", 0.1));
    D.UI.push_back(ParamUI("ForceGrav_______", -9.81));
    D.UI.push_back(ParamUI("StiffBox________", 1.e3));
    D.UI.push_back(ParamUI("StiffCollision__", 1.e1));
    D.UI.push_back(ParamUI("StiffEdgeLength_", 1.e1));
    D.UI.push_back(ParamUI("StiffTriArea____", 1.e3));
    D.UI.push_back(ParamUI("StiffTetVolume__", 1.e6));
    D.UI.push_back(ParamUI("StiffMouseBall__", 1.e2));
    D.UI.push_back(ParamUI("RandPerturb_____", 1.e-6));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("ColorMode_______", 2));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("VisuSimple______", -1));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("TestParamPBD_0__", 0.0));
    D.UI.push_back(ParamUI("TestParamPBD_1__", 0.0));
    D.UI.push_back(ParamUI("TestParamPBD_2__", 0.0));
    D.UI.push_back(ParamUI("TestParamPBD_3__", 0.0));
    D.UI.push_back(ParamUI("TestParamPBD_4__", 0.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));

    D.displayModeLabel[1]= "Nodes";
    D.displayModeLabel[2]= "Edges";
    D.displayModeLabel[3]= "Tris";
    D.displayModeLabel[4]= "Tets";
    D.displayModeLabel[5]= "Normals";
    D.displayModeLabel[6]= "Mouse";
    D.displayMode[2]= false;
    D.displayMode[4]= false;
    D.displayMode[5]= false;
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool PosiBasedDynam::CheckAlloc() {
  if (D.UI[NumParticl______].hasChanged()) isAllocated= false;
  if (D.UI[AddTetMeshID____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool PosiBasedDynam::CheckRefresh() {
  if (D.UI[RadParticl______].hasChanged()) isRefreshed= false;
  if (D.UI[MaterialDensity_].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void PosiBasedDynam::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get domain dimensions
  D.boxMin= {0.5f - 0.5f * D.UI[DomainX_________].F(), 0.5f - 0.5f * D.UI[DomainY_________].F(), 0.5f - 0.5f * D.UI[DomainZ_________].F()};
  D.boxMax= {0.5f + 0.5f * D.UI[DomainX_________].F(), 0.5f + 0.5f * D.UI[DomainY_________].F(), 0.5f + 0.5f * D.UI[DomainZ_________].F()};

  // Clear data
  Pos.clear();
  Edge.clear();
  Tri.clear();
  Tet.clear();

  // Add random individual nodes
  std::default_random_engine rng(0);
  std::uniform_real_distribution<float> distribX(D.boxMin[0], D.boxMax[0]), distribY(D.boxMin[1], D.boxMax[1]), distribZ(D.boxMin[2], D.boxMax[2]);
  for (int k= 0; k < D.UI[NumParticl______].I(); k++) {
    Pos.push_back({distribX(rng), distribY(rng), distribZ(rng)});
  }

  // Add preset tet meshes
  if (D.UI[AddTetMeshID____].I() == 0) {
    const int offset= (int)Pos.size();
    Pos.push_back({0.3f, 0.3f, 0.3f});
    Pos.push_back({0.7f, 0.3f, 0.3f});
    Pos.push_back({0.3f, 0.7f, 0.3f});
    Pos.push_back({0.3f, 0.3f, 0.7f});
    Edge.push_back({offset+0, offset+1});
    Edge.push_back({offset+0, offset+2});
    Edge.push_back({offset+0, offset+3});
    Edge.push_back({offset+1, offset+2});
    Edge.push_back({offset+1, offset+3});
    Edge.push_back({offset+2, offset+3});
    Tri.push_back({offset+0, offset+2, offset+1});
    Tri.push_back({offset+0, offset+1, offset+3});
    Tri.push_back({offset+0, offset+3, offset+2});
    Tri.push_back({offset+1, offset+2, offset+3});
    Tet.push_back({offset+0, offset+1, offset+2, offset+3});
  }
  else if (D.UI[AddTetMeshID____].I() > 0) {
    const int offset= (int)Pos.size();
    std::vector<std::array<float, 3>> loadedPos;
    std::vector<std::array<int, 4>> loadedTet;
    if (D.UI[AddTetMeshID____].I() == 1) FileInput::LoadTetMeshMSHFile("./FileInput/TetMesh/Bunny.msh", loadedPos, loadedTet, true);
    if (D.UI[AddTetMeshID____].I() == 2) FileInput::LoadTetMeshMSHFile("./FileInput/TetMesh/Kitten.msh", loadedPos, loadedTet, true);
    if (D.UI[AddTetMeshID____].I() == 3) FileInput::LoadTetMeshMSHFile("./FileInput/TetMesh/AssEars.msh", loadedPos, loadedTet, true);
    if (!loadedPos.empty() && !loadedTet.empty()) {
      // Add the loaded nodes
      for (unsigned int k= 0; k < loadedPos.size(); k++) {
        Pos.push_back(loadedPos[k]);
      }
      // Add the loaded tetrahedra
      for (unsigned int k= 0; k < loadedTet.size(); k++) {
        Tet.push_back({offset + loadedTet[k][0], offset + loadedTet[k][1], offset + loadedTet[k][2], offset + loadedTet[k][3]});
      }
      // Create and add the bars
      for (unsigned int k= 0; k < loadedTet.size(); k++) {
        if (loadedTet[k][0] < loadedTet[k][1]) Edge.push_back({offset + loadedTet[k][0], offset + loadedTet[k][1]});
        if (loadedTet[k][0] < loadedTet[k][2]) Edge.push_back({offset + loadedTet[k][0], offset + loadedTet[k][2]});
        if (loadedTet[k][0] < loadedTet[k][3]) Edge.push_back({offset + loadedTet[k][0], offset + loadedTet[k][3]});
        if (loadedTet[k][1] < loadedTet[k][2]) Edge.push_back({offset + loadedTet[k][1], offset + loadedTet[k][2]});
        if (loadedTet[k][1] < loadedTet[k][3]) Edge.push_back({offset + loadedTet[k][1], offset + loadedTet[k][3]});
        if (loadedTet[k][2] < loadedTet[k][3]) Edge.push_back({offset + loadedTet[k][2], offset + loadedTet[k][3]});
      }
      // Build the list of all faces in the tetrahedral mesh
      std::vector<std::array<int, 3>> TmpTri;
      for (unsigned int k= 0; k < loadedTet.size(); k++) {
        TmpTri.push_back({offset + loadedTet[k][0], offset + loadedTet[k][2], offset + loadedTet[k][1]});
        TmpTri.push_back({offset + loadedTet[k][0], offset + loadedTet[k][1], offset + loadedTet[k][3]});
        TmpTri.push_back({offset + loadedTet[k][0], offset + loadedTet[k][3], offset + loadedTet[k][2]});
        TmpTri.push_back({offset + loadedTet[k][1], offset + loadedTet[k][2], offset + loadedTet[k][3]});
      }
      // Filter the list to only insert the external faces of the tetrahedral mesh
      std::map<std::array<int, 3>, char> encounters;
      for (unsigned int k= 0; k < TmpTri.size(); k++) {
        std::array<int, 3> curTri= TmpTri[k];
        std::sort(curTri.begin(), curTri.end());
        encounters[curTri]++;
      }
      for (unsigned int k= 0; k < TmpTri.size(); k++) {
        std::array<int, 3> curTri= TmpTri[k];
        std::sort(curTri.begin(), curTri.end());
        if (encounters[curTri] == 1) Tri.push_back(TmpTri[k]);
      }
    }
  }

  // Set the previous position
  PosOld= Pos;

  // Set the initial velocity
  Vel= std::vector<Vec::Vec3<float>>(Pos.size(), {0.0f, 0.0f, 0.0f});
}


// Refresh the project
void PosiBasedDynam::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Compute the rest length for bars
  EdgeLen.clear();
  for (int k= 0; k < (int)Edge.size(); k++)
    EdgeLen.push_back(GetEdgeLength(Pos[Edge[k][0]], Pos[Edge[k][1]]));

  // Compute the rest area for triangles
  TriArea.clear();
  for (int k= 0; k < (int)Tri.size(); k++)
    TriArea.push_back(GetTriArea(Pos[Tri[k][0]], Pos[Tri[k][1]], Pos[Tri[k][2]]));

  // Compute the rest volumes for tetrahedra
  TetVol.clear();
  for (int k= 0; k < (int)Tet.size(); k++)
    TetVol.push_back(GetTetVolume(Pos[Tet[k][0]], Pos[Tet[k][1]], Pos[Tet[k][2]], Pos[Tet[k][3]]));
  
  // Compute masses and inverses
  Mass= std::vector<float>(Pos.size(), 0.0f);
  for (int k= 0; k < D.UI[NumParticl______].I(); k++)
    Mass[k]= D.UI[MaterialDensity_].F() * (4.0f/3.0f) * std::numbers::pi * std::pow(D.UI[RadParticl______].F(), 3.0f);
  for (int k= 0; k < (int)Tet.size(); k++)
    for (int idx : Tet[k])
      Mass[idx]+= 0.25f * D.UI[MaterialDensity_].F() * TetVol[k];

  MassInv.clear();
  for (int k= 0; k < (int)Pos.size(); k++)
    MassInv.push_back(1.0f / Mass[k]);
}


// Handle UI parameter change
void PosiBasedDynam::ParamChange() {
  if (D.UI[DomainX_________].hasChanged() || D.UI[DomainY_________].hasChanged() || D.UI[DomainZ_________].hasChanged()) {
    D.boxMin= {0.5f - 0.5f * D.UI[DomainX_________].F(), 0.5f - 0.5f * D.UI[DomainY_________].F(), 0.5f - 0.5f * D.UI[DomainZ_________].F()};
    D.boxMax= {0.5f + 0.5f * D.UI[DomainX_________].F(), 0.5f + 0.5f * D.UI[DomainY_________].F(), 0.5f + 0.5f * D.UI[DomainZ_________].F()};
  }
}


// Handle keypress
void PosiBasedDynam::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  // Stability check by entirely collapsing the entire geometry to a single immobile point
  if (D.keyLetterUpperCase == 'C') {
    Pos= std::vector<Vec::Vec3<float>>(Pos.size(), {0.5f, 0.5f, 0.5f});
    PosOld= std::vector<Vec::Vec3<float>>(Pos.size(), {0.5f, 0.5f, 0.5f});
    Vel= std::vector<Vec::Vec3<float>>(Pos.size(), {0.0f, 0.0f, 0.0f});
  }
}


// Handle mouse action
void PosiBasedDynam::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void PosiBasedDynam::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  TimeIntegrate();
}


// Draw the project
void PosiBasedDynam::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;
  
  DrawScene();
}
