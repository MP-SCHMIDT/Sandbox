#include "NBodyGravDynam.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <format>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Type/Vec.hpp"
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
NBodyGravDynam::NBodyGravDynam() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void NBodyGravDynam::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("BodyCount_______", 100));
    D.UI.push_back(ParamUI("BodyStartLayout_", 0));
    D.UI.push_back(ParamUI("BodyInitVel_____", 0.1));
    D.UI.push_back(ParamUI("BodyRadius______", 0.01));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("Domain2D________", 0));
    D.UI.push_back(ParamUI("DomainTaurusPos_", 0));
    D.UI.push_back(ParamUI("DomainTaurusFor_", 0));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("TreeMaxDepth____", 8));
    D.UI.push_back(ParamUI("TreeTolRatio____", 0.5));
    D.UI.push_back(ParamUI("TreeShowMin_____", 0));
    D.UI.push_back(ParamUI("TreeShowMax_____", 10));
    D.UI.push_back(ParamUI("TreeShowEmpty___", 0));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("SimuStepMode____", 1));
    D.UI.push_back(ParamUI("SimuForceMode___", 1));
    D.UI.push_back(ParamUI("SimuTimeStep____", 0.005));
    D.UI.push_back(ParamUI("SimuTotGravity__", 1.0));
    D.UI.push_back(ParamUI("______________03", NAN));
    D.UI.push_back(ParamUI("ColorMode_______", 2));
    D.UI.push_back(ParamUI("ColorFactor_____", 0.1));
    D.UI.push_back(ParamUI("______________04", NAN));
    D.UI.push_back(ParamUI("TestParamNBS_00_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_01_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_02_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_03_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_04_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_05_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_06_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_07_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_08_", 0.0));
    D.UI.push_back(ParamUI("TestParamNBS_09_", 0.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));

    D.displayModeLabel[1]= "Bodies";
    D.displayModeLabel[2]= "BodiesSpheres";
    D.displayModeLabel[3]= "Octree";
    D.displayModeLabel[4]= "OctreeMass";
    #ifdef TESTING_DISPLAY_FORCES_VECTORS
    D.displayModeLabel[5]= "ApproxForce";
    D.displayModeLabel[6]= "ApproxSource";
    #endif
    D.displayMode[2]= false;
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool NBodyGravDynam::CheckAlloc() {
  if (D.UI[BodyCount_______].hasChanged()) isAllocated= false;
  if (D.UI[BodyStartLayout_].hasChanged()) isAllocated= false;
  if (D.UI[BodyInitVel_____].hasChanged()) isAllocated= false;
  if (D.UI[Domain2D________].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool NBodyGravDynam::CheckRefresh() {
  if (D.UI[TreeMaxDepth____].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void NBodyGravDynam::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Set the box domain
  D.boxMin= {0.0f, 0.0f, 0.0f};
  D.boxMax= {1.0f, 1.0f, 1.0f};

  // Get UI parameters
  N= std::max(D.UI[BodyCount_______].I(), 1);

  // Allocate data
  Pos= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  Vel= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  Acc= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  For= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));

  // Initialize bodies
  if (D.UI[BodyStartLayout_].I() == 0) {
    // Random positions and velocities
    for (unsigned int k0= 0; k0 < N; k0++) {
      Pos[k0].set(Random::Val(0.0f, 1.0f), Random::Val(0.0f, 1.0f), Random::Val(0.0f, 1.0f));
      Vel[k0]= Vec::Vec3<float>(Random::Val(-1.0f, 1.0f), Random::Val(-1.0f, 1.0f), Random::Val(-1.0f, 1.0f)) * D.UI[BodyInitVel_____].F();
    }
  }
  else if (D.UI[BodyStartLayout_].I() == 1) {
    // Sphere with angular velocity
    const Vec::Vec3<float> center(0.5f, 0.5f, 0.5f);
    const Vec::Vec3<float> spinAxis(1.0f, 0.0f, 0.0f);
    for (unsigned int k0= 0; k0 < N; k0++) {
      for (int idxDart= 0; idxDart < 100; idxDart++) {
        Pos[k0].set(Random::Val(0.25f, 0.75f), Random::Val(0.25f, 0.75f), Random::Val(0.25f, 0.75f));
        if ((center - Pos[k0]).norm() < 0.25f) break;
      }
      if ((center - Pos[k0]).norm() > 0.0f) {
        const Vec::Vec3<float> tangent= spinAxis.cross(Pos[k0]-center) / spinAxis.cross(Pos[k0]-center).norm();
        Vel[k0]= tangent * D.UI[BodyInitVel_____].F() * Random::Val(0.5f, 1.5f) * (center - Pos[k0]).norm();
      } 
    }
  }

  // Ensure positions are in the valid domain
  ApplyPosBC();
}


// Refresh the project
void NBodyGravDynam::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Initialize the tree for visualization purposes
  BuildTree();
}


// Handle UI parameter change
void NBodyGravDynam::ParamChange() {
}


// Handle keypress
void NBodyGravDynam::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Handle mouse action
void NBodyGravDynam::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void NBodyGravDynam::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Get UI parameters
  const float dt= D.UI[SimuTimeStep____].F();

  // Explicit Euler integration
  if (D.UI[SimuStepMode____].I() == 0) {
    BuildTree();
    ComputeForces();  // f(x₀)
    for (unsigned int k0= 0; k0 < N; k0++) {
      Acc[k0]= For[k0];                 // a₁ = f(x₀) / m
      Vel[k0]= Vel[k0] + dt * Acc[k0];  // v₁ = v₀ + Δt a₁
      Pos[k0]= Pos[k0] + dt * Vel[k0];  // x₁ = x₀ + Δt v₁
    }
    ApplyPosBC();
  }

  // Explicit Velocity Verlet integration
  if (D.UI[SimuStepMode____].I() == 1) {
    for (unsigned int k0= 0; k0 < N; k0++) {
      Pos[k0]= Pos[k0] + Vel[k0] * dt + 0.5f * Acc[k0] * dt * dt;  // x₁ = x₀ + Δt v₀ + 0.5 * a₀ * Δt²
    }
    ApplyPosBC();
    BuildTree();
    ComputeForces();  // f(x₁)
    for (unsigned int k0= 0; k0 < N; k0++) {
      Vec::Vec3<float> AccOld= Acc[k0];
      Acc[k0]= For[k0];                                   // a₁ = f(x₁) / m
      Vel[k0]= Vel[k0] + 0.5f * (AccOld + Acc[k0]) * dt;  // v₁ = v₀ + 0.5 * (a₀ + a₁) * Δt
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("TimeAnim %f ", Timer::PopTimer());
  if (D.UI[VerboseLevel____].I() >= 1) printf("\n");
}


// Draw the project
void NBodyGravDynam::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the particles
  if (D.displayMode[1] || D.displayMode[2]) {
    if (D.displayMode[1]) {
      glPointSize(1000.0f * D.UI[BodyRadius______].F());
      glEnable(GL_POINT_SMOOTH);
      glBegin(GL_POINTS);
    }
    else if (D.displayMode[2]) {
      glEnable(GL_LIGHTING);
    }
    for (unsigned int k0= 0; k0 < N; k0++) {
      float r, g, b;
      if (D.UI[ColorMode_______].I() == 0) {
        Colormap::RatioToJetBrightSmooth(Vel[k0].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      }
      if (D.UI[ColorMode_______].I() == 1) {
        r= 0.5f + Vel[k0][0] * D.UI[ColorFactor_____].F();
        g= 0.5f + Vel[k0][1] * D.UI[ColorFactor_____].F();
        b= 0.5f + Vel[k0][2] * D.UI[ColorFactor_____].F();
      }
      if (D.UI[ColorMode_______].I() == 2) {
        Colormap::RatioToPlasma(For[k0].norm() * D.UI[ColorFactor_____].F(), r, g, b);
      }
      glColor3f(r, g, b);

      if (D.displayMode[1]) {
        glVertex3fv(Pos[k0].array());
      }
      else if (D.displayMode[2]) {
        glPushMatrix();
        glTranslatef(Pos[k0][0], Pos[k0][1], Pos[k0][2]);
        glScalef(D.UI[BodyRadius______].F(), D.UI[BodyRadius______].F(), D.UI[BodyRadius______].F());
        glutSolidSphere(1.0, 32, 16);
        glPopMatrix();
      }
    }
    if (D.displayMode[1]) {
      glEnd();
      glDisable(GL_POINT_SMOOTH);
      glPointSize(1.0f);
    }
    else if (D.displayMode[2]) {
      glDisable(GL_LIGHTING);
    }
  }

  // Draw the octree
  if (D.displayMode[3] || D.displayMode[4]) {
    for (NBodyGravDynam::OctreeNode Cell : Tree) {
      if (Cell.Depth >= D.UI[TreeShowMin_____].I() && Cell.Depth <= D.UI[TreeShowMax_____].I()) {
        if (D.UI[TreeShowEmpty___].I() || Cell.Count > 0) {
          if (D.displayMode[3]) {
            float r, g, b;
            Colormap::RatioToRainbow((float)Cell.Depth * 0.1f, r, g, b);
            glColor3f(r, g, b);
            glPushMatrix();
            glTranslatef(Cell.Center[0], Cell.Center[1], Cell.Center[2]);
            glutWireCube(Cell.Size);
            glPopMatrix();
          }
          if (D.displayMode[4]) {
            float r, g, b;
            Colormap::RatioToTurbo((float)Cell.Count * D.UI[ColorFactor_____].F(), r, g, b);
            glColor3f(r, g, b);
            glPushMatrix();
            glTranslatef(Cell.CenterMass[0], Cell.CenterMass[1], Cell.CenterMass[2]);
            glutWireCube(Cell.Size * 0.2f);
            glPopMatrix();
          }
        }
      }
    }
  }

  // Draw the approximation of forces for the chosen particle
  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  if (D.displayMode[5]) {
    if (D.UI[TestParamNBS_00_].I() >= 0 && D.UI[TestParamNBS_00_].I() < (int)N) {
      glLineWidth(2.0f);
      glBegin(GL_LINES);
      Vec::Vec3<float> p0= Pos[D.UI[TestParamNBS_00_].I()];
      for (unsigned int k= 0; k < ContribPos.size(); k++) {
        float r, g, b;
        Colormap::RatioToTurbo((float)ContribCount[k] * D.UI[ColorFactor_____].F(), r, g, b);
        glColor3f(r, g, b);
        glVertex3fv(p0.array());
        glVertex3fv(ContribPos[k].array());
      }
      glEnd();
      glLineWidth(1.0f);
    }
  }
  #endif

  
  // Draw the contributing cells
  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  if (D.displayMode[6]) {
    for (unsigned int idxCell : ContribCell) {
      const NBodyGravDynam::OctreeNode Cell= Tree[idxCell];
      if (Cell.Depth >= D.UI[TreeShowMin_____].I() && Cell.Depth <= D.UI[TreeShowMax_____].I()) {
        float r, g, b;
        Colormap::RatioToRainbow((float)Cell.Depth * 0.1f, r, g, b);
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslatef(Cell.Center[0], Cell.Center[1], Cell.Center[2]);
        glutWireCube(Cell.Size);
        glPopMatrix();
      }
    }
  }
  #endif

  D.Status.clear();
  D.Status.resize(1);
  D.Status[0]= std::format("Tree cells:{}", (int)Tree.size());
}


void NBodyGravDynam::AddTreeNode(const float iCenterX,
                                 const float iCenterY,
                                 const float iCenterZ,
                                 const float iSize,
                                 const unsigned char iDepth) {
  const unsigned int idxCell= Tree.size();
  Tree.push_back(NBodyGravDynam::OctreeNode());
  Tree[idxCell].Center.set(iCenterX, iCenterY, iCenterZ);
  Tree[idxCell].Size= iSize;
  Tree[idxCell].Depth= iDepth;
  Tree[idxCell].CenterMass.set(0.0f, 0.0f, 0.0f);
  Tree[idxCell].Count= 0;
  Tree[idxCell].Child= 0;
}


void NBodyGravDynam::BuildTree() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Initialize the tree with the root cell
  Tree.clear();
  AddTreeNode(0.5f, 0.5f, 0.5f, 1.0f, 0);

  // Sequentially add the particles in the tree, adding cells as required
  for (unsigned int k0= 0; k0 < N; k0++) {
    // Add the particle to the root cell mass and count
    Tree[0].CenterMass= Tree[0].CenterMass + Pos[k0];
    Tree[0].Count++;
    
    // Start the tree search
    unsigned int idxCell= 0;
    while (Tree[idxCell].Depth < D.UI[TreeMaxDepth____].I()) {
      // Add child cells if not already present
      if (Tree[idxCell].Child == 0) {
        Tree[idxCell].Child= Tree.size();
        const float childSize= 0.5f * Tree[idxCell].Size;
        const float centerXN= Tree[idxCell].Center[0] - 0.5f * childSize;
        const float centerXP= Tree[idxCell].Center[0] + 0.5f * childSize;
        const float centerYN= Tree[idxCell].Center[1] - 0.5f * childSize;
        const float centerYP= Tree[idxCell].Center[1] + 0.5f * childSize;
        const float centerZN= Tree[idxCell].Center[2] - 0.5f * childSize;
        const float centerZP= Tree[idxCell].Center[2] + 0.5f * childSize;
        const unsigned char childDepth= Tree[idxCell].Depth + 1;
        AddTreeNode(centerXN, centerYN, centerZN, childSize, childDepth);
        AddTreeNode(centerXN, centerYN, centerZP, childSize, childDepth);
        AddTreeNode(centerXN, centerYP, centerZN, childSize, childDepth);
        AddTreeNode(centerXN, centerYP, centerZP, childSize, childDepth);
        AddTreeNode(centerXP, centerYN, centerZN, childSize, childDepth);
        AddTreeNode(centerXP, centerYN, centerZP, childSize, childDepth);
        AddTreeNode(centerXP, centerYP, centerZN, childSize, childDepth);
        AddTreeNode(centerXP, centerYP, centerZP, childSize, childDepth);
      }

      // Find child cell where the body belongs
      unsigned int idxChild= Tree[idxCell].Child;
      if (Pos[k0][0] > Tree[idxCell].Center[0]) idxChild+= 4;
      if (Pos[k0][1] > Tree[idxCell].Center[1]) idxChild+= 2;
      if (Pos[k0][2] > Tree[idxCell].Center[2]) idxChild+= 1;

      // Update the mass and count of the child cell
      Tree[idxChild].CenterMass= Tree[idxChild].CenterMass + Pos[k0];
      Tree[idxChild].Count++;

      // Set the child cell as new cell for the next depth
      idxCell= idxChild;
    }
  }

  // Compute average center of mass for each tree cell
  for (unsigned int idxCell= 0; idxCell < Tree.size(); idxCell++)
    if (Tree[idxCell].Count > 0)
      Tree[idxCell].CenterMass= Tree[idxCell].CenterMass/float(Tree[idxCell].Count);

  if (D.UI[VerboseLevel____].I() >= 1) printf("TimeTreeBuild %f ", Timer::PopTimer());
}


void NBodyGravDynam::ComputeForces() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  const float gravity= D.UI[SimuTotGravity__].F() / float(D.UI[BodyCount_______].I());
  const float minRadSqr= (2.0f * D.UI[BodyRadius______].F()) * (2.0f * D.UI[BodyRadius______].F());

  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  if (D.displayMode[5]) {
    ContribPos.clear();
    ContribCount.clear();
    ContribCell.clear();
  }
  #endif

  if (D.UI[SimuForceMode___].I() == 0) {
    #pragma omp parallel for
    for (unsigned int k0= 0; k0 < N; k0++) {
      const Vec::Vec3<float> p0(Pos[k0]);
      For[k0].set(0.0f, 0.0f, 0.0f);
      for (unsigned int k1= 0; k1 < N; k1++) {
        if (k0 == k1) continue;
  
        Vec::Vec3<float> p1(Pos[k1]);
        if (D.UI[DomainTaurusFor_].I() == 1) {
          if      ((p1[0] - p0[0]) >  0.5f) p1[0]-= 1.0f;
          else if ((p1[0] - p0[0]) < -0.5f) p1[0]+= 1.0f;
          if      ((p1[1] - p0[1]) >  0.5f) p1[1]-= 1.0f;
          else if ((p1[1] - p0[1]) < -0.5f) p1[1]+= 1.0f;
          if      ((p1[2] - p0[2]) >  0.5f) p1[2]-= 1.0f;
          else if ((p1[2] - p0[2]) < -0.5f) p1[2]+= 1.0f;
        }
  
        const Vec::Vec3<float> vec(p1 - p0);
        For[k0]+= (vec/vec.norm()) * gravity / std::max(vec.normSquared(), minRadSqr);
        
        #ifdef TESTING_DISPLAY_FORCES_VECTORS
        if (D.displayMode[5] && (int)k0 == D.UI[TestParamNBS_00_].I()) {
          ContribPos.push_back(p1);
          ContribCount.push_back(1);
        }
        #endif
      }
    }
  }


  if (D.UI[SimuForceMode___].I() == 1) {
    const float treeTolSqr=  D.UI[TreeTolRatio____].F() * D.UI[TreeTolRatio____].F();
    #pragma omp parallel for
    for (unsigned int k0= 0; k0 < N; k0++) {
      const Vec::Vec3<float> p0(Pos[k0]);
      For[k0].set(0.0f, 0.0f, 0.0f);
  
      std::vector<unsigned int> Q(1, 0);
      while (Q.size() > 0) {
        #ifdef TESTING_DISPLAY_FORCES_VECTORS
        const unsigned int idxCell= Q[Q.size()-1];
        #endif
        const NBodyGravDynam::OctreeNode cell= Tree[Q[Q.size()-1]];
        Q.pop_back();
        
        Vec::Vec3<float> p1(cell.CenterMass);
        if (D.UI[DomainTaurusFor_].I() == 1) {
          if      ((p1[0] - p0[0]) >  0.5f) p1[0]-= 1.0f;
          else if ((p1[0] - p0[0]) < -0.5f) p1[0]+= 1.0f;
          if      ((p1[1] - p0[1]) >  0.5f) p1[1]-= 1.0f;
          else if ((p1[1] - p0[1]) < -0.5f) p1[1]+= 1.0f;
          if      ((p1[2] - p0[2]) >  0.5f) p1[2]-= 1.0f;
          else if ((p1[2] - p0[2]) < -0.5f) p1[2]+= 1.0f;
        }

        const Vec::Vec3<float> vec(p1 - p0);
        if (cell.Child == 0 || (cell.Size*cell.Size) / vec.normSquared() < treeTolSqr) {
          if (vec[0] != 0.0f || vec[2] != 0.0f || vec[2] != 0.0f) {
            For[k0]+= (vec/vec.norm()) * gravity * float(cell.Count) / std::max(vec.normSquared(), minRadSqr);
            #ifdef TESTING_DISPLAY_FORCES_VECTORS
            if (D.displayMode[5] && (int)k0 == D.UI[TestParamNBS_00_].I()) {
              ContribPos.push_back(p1);
              ContribCount.push_back(cell.Count);
              ContribCell.push_back(idxCell);
            }
            #endif
          }
        }
        else {
          for(unsigned int idxChild= cell.Child; idxChild < cell.Child+8; idxChild++) {
            if (Tree[idxChild].Count > 0) {
              Q.push_back(idxChild);
            }
          }
        }
      }
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("TimeCompForces %f ", Timer::PopTimer());
}


void NBodyGravDynam::ApplyPosBC() {
  if (D.UI[Domain2D________].I() == 1) {
    for (unsigned int k0= 0; k0 < N; k0++)
      Acc[k0][0]= Vel[k0][0]= Pos[k0][0]= 0.51f;
  }

  if (D.UI[DomainTaurusPos_].I() == 1) {
    for (unsigned int k0= 0; k0 < N; k0++) {
      if (Pos[k0][0] < 0.0f || Pos[k0][0] > 1.0f) Pos[k0][0]= Pos[k0][0] - std::floor(Pos[k0][0]);
      if (Pos[k0][1] < 0.0f || Pos[k0][1] > 1.0f) Pos[k0][1]= Pos[k0][1] - std::floor(Pos[k0][1]);
      if (Pos[k0][2] < 0.0f || Pos[k0][2] > 1.0f) Pos[k0][2]= Pos[k0][2] - std::floor(Pos[k0][2]);
    }
  }
}
