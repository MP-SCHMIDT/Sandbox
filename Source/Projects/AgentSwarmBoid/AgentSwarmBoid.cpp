#include "AgentSwarmBoid.hpp"


// Standard lib
#include <vector>

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
AgentSwarmBoid::AgentSwarmBoid() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void AgentSwarmBoid::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("Constrain2D_____", 0));
    D.UI.push_back(ParamUI("PopSize_________", 300));
    D.UI.push_back(ParamUI("PopTypes________", 3));
    D.UI.push_back(ParamUI("TimeStep________", 0.05));
    D.UI.push_back(ParamUI("SizeView________", 0.15));
    D.UI.push_back(ParamUI("SizeBody________", 0.05));
    D.UI.push_back(ParamUI("CoeffSep________", 0.15));
    D.UI.push_back(ParamUI("CoeffAli________", 0.05));
    D.UI.push_back(ParamUI("CoeffCoh________", 0.10));
    D.UI.push_back(ParamUI("CoeffEat________", 0.10));
    D.UI.push_back(ParamUI("CoeffRun________", 0.15));
    D.UI.push_back(ParamUI("CoeffOri________", 0.02));
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
bool AgentSwarmBoid::CheckAlloc() {
  if (D.UI[PopSize_________].hasChanged()) isAllocated= false;
  if (D.UI[PopTypes________].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool AgentSwarmBoid::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void AgentSwarmBoid::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  NbAgents= std::max(D.UI[PopSize_________].I(), 1);
  NbTypes= std::max(D.UI[PopTypes________].I(), 1);

  // Allocate data
  Pos= std::vector<Vec::Vec3<float>>(NbAgents);
  Vel= std::vector<Vec::Vec3<float>>(NbAgents);
  Typ= std::vector<int>(NbAgents);
}


// Refresh the project
void AgentSwarmBoid::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Initialize population with random positions, velocities and types
  for (int k= 0; k < NbAgents; k++) {
    Pos[k].set(Random::Val(0.0f, 1.0f), Random::Val(0.0f, 1.0f), Random::Val(0.0f, 1.0f));
    Vel[k].set(Random::Val(-0.2f, 0.2f), Random::Val(-0.2f, 0.2f), Random::Val(-0.2f, 0.2f));
    Typ[k]= Random::Val(0, NbTypes - 1);
  }

  // Optionally constrain to 2D
  if (D.UI[Constrain2D_____].B()) {
    for (int k0= 0; k0 < NbAgents; k0++) {
      Pos[k0][0]= 0.5;
      Vel[k0][0]= 0.0;
    }
  }
}


// Handle keypress
void AgentSwarmBoid::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Handle mouse action
void AgentSwarmBoid::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
}


// Animate the project
void AgentSwarmBoid::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Optionally constrain to 2D
  if (D.UI[Constrain2D_____].B()) {
    for (int k0= 0; k0 < NbAgents; k0++) {
      Pos[k0][0]= 0.5;
      Vel[k0][0]= 0.0;
    }
  }

  // Compute the forces
  std::vector<Vec::Vec3<float>> velocityChange(NbAgents);
#pragma omp parallel for
  for (int k0= 0; k0 < NbAgents; k0++) {
    int countSep= 0;
    int countAli= 0;
    int countEat= 0;
    int countRun= 0;
    Vec::Vec3<float> sep, ali, coh, eat, run, ori;
    for (int k1= 0; k1 < NbAgents; k1++) {
      if ((Pos[k0] - Pos[k1]).normSquared() < D.UI[SizeView________].F() * D.UI[SizeView________].F()) {
        if (k0 != k1) {
          if ((Pos[k0] - Pos[k1]).normSquared() < D.UI[SizeBody________].F() * D.UI[SizeBody________].F()) {
            sep+= Pos[k0] - Pos[k1];
            countSep++;
          }
          if (Typ[k0] == Typ[k1]) {
            ali+= Vel[k1];
            coh= Pos[k1] - Pos[k0];
            countAli++;
          }
          if (Typ[k0] == (Typ[k1] + NbTypes - 1) % NbTypes) {
            eat= Pos[k1] - Pos[k0];
            countEat++;
          }
          if (Typ[k0] == (Typ[k1] + 1) % NbTypes) {
            run= Pos[k0] - Pos[k1];
            countRun++;
          }
        }
      }
    }
    if (countSep > 0) sep/= countSep;
    if (countAli > 0) ali/= countAli;
    if (countAli > 0) coh/= countAli;
    if (countEat > 0) eat/= countEat;
    if (countRun > 0) run/= countRun;

    Vec::Vec3<float> center(0.5f, 0.5f, 0.5f);
    ori= center - Pos[k0];

    velocityChange[k0].set(0.0f, 0.0f, 0.0f);
    velocityChange[k0]+= D.UI[CoeffSep________].F() * sep;
    velocityChange[k0]+= D.UI[CoeffAli________].F() * ali;
    velocityChange[k0]+= D.UI[CoeffCoh________].F() * coh;
    velocityChange[k0]+= D.UI[CoeffEat________].F() * eat;
    velocityChange[k0]+= D.UI[CoeffRun________].F() * run;
    velocityChange[k0]+= D.UI[CoeffOri________].F() * ori;
  }

  // Apply the forces
  for (int k= 0; k < NbAgents; k++) {
    Vel[k]+= velocityChange[k];
    float l= Vel[k].norm();
    if (l > 0.3f) {
      Vel[k]= 0.3f * Vel[k] / l;
    }
    Pos[k]= Pos[k] + Vel[k] * D.UI[TimeStep________].F();
  }
}


// Draw the project
void AgentSwarmBoid::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the agents
  if (D.displayMode1) {
    for (int k= 0; k < NbAgents; k++) {
      Vec::Vec3<float> front= Vel[k].normalized();
      Vec::Vec3<float> u(1.0f, 0.0f, 0.0f);
      if (u.dot(front) == 0.0f) u.set(0.0f, 1.0f, 0.0f);
      Vec::Vec3<float> left= front.cross(u).normalized();
      Vec::Vec3<float> up= front.cross(left).normalized();

      Vec::Vec3<float> p1= Pos[k] + 0.05f * 2.0f * (+0.5f * front + 0.00f * left + 0.0f * up);
      Vec::Vec3<float> p2= Pos[k] + 0.05f * 2.0f * (+0.0f * front - 0.07f * left + 0.1f * up);
      Vec::Vec3<float> p3= Pos[k] + 0.05f * 2.0f * (+0.0f * front + 0.07f * left + 0.1f * up);
      Vec::Vec3<float> p4= Pos[k] + 0.05f * 2.0f * (+0.0f * front + 0.00f * left - 0.1f * up);
      Vec::Vec3<float> p5= Pos[k] + 0.05f * 2.0f * (-0.4f * front + 0.00f * left + 0.0f * up);
      Vec::Vec3<float> p6= Pos[k] + 0.05f * 2.0f * (-0.6f * front + 0.00f * left - 0.2f * up);
      Vec::Vec3<float> p7= Pos[k] + 0.05f * 2.0f * (-0.6f * front + 0.00f * left + 0.2f * up);

      float r= 0.5f, g= 0.5f, b= 0.5f;
      if (NbTypes > 1)
        Colormap::RatioToRainbow((float)Typ[k] / (float)(NbTypes - 1), r, g, b);

      glBegin(GL_TRIANGLES);
      glColor3f(r, g, b);
      glVertex3fv(p1.array());
      glVertex3fv(p2.array());
      glVertex3fv(p3.array());
      glVertex3fv(p1.array());
      glVertex3fv(p4.array());
      glVertex3fv(p2.array());
      glVertex3fv(p1.array());
      glVertex3fv(p3.array());
      glVertex3fv(p4.array());
      glVertex3fv(p2.array());
      glVertex3fv(p4.array());
      glVertex3fv(p5.array());
      glVertex3fv(p3.array());
      glVertex3fv(p2.array());
      glVertex3fv(p5.array());
      glVertex3fv(p3.array());
      glVertex3fv(p5.array());
      glVertex3fv(p4.array());
      glVertex3fv(p5.array());
      glVertex3fv(p6.array());
      glVertex3fv(p7.array());
      glEnd();
    }
  }

  // Draw the agents velocity vectors
  if (D.displayMode2) {
    glBegin(GL_LINES);
    for (int k= 0; k < NbAgents; k++) {
      glColor3f(0.0f, 0.0f, 1.0f);
      glVertex3fv(Pos[k].array());
      glColor3f(1.0f, 0.0f, 0.0f);
      glVertex3fv((Pos[k] + Vel[k]).array());
    }
    glEnd();
  }
}
