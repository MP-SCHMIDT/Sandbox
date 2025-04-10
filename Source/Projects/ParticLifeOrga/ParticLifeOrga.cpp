#include "ParticLifeOrga.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Type/Field.hpp"
#include "Type/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
ParticLifeOrga::ParticLifeOrga() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void ParticLifeOrga::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("NbTypes_________", 6));
    D.UI.push_back(ParamUI("NbPartic________", 800));
    D.UI.push_back(ParamUI("TimeStep________", 0.05));
    D.UI.push_back(ParamUI("NumIntegMode____", 0));
    D.UI.push_back(ParamUI("Use2DProj_______", 1));
    D.UI.push_back(ParamUI("DomainX_________", 1.0));
    D.UI.push_back(ParamUI("DomainY_________", 1.0));
    D.UI.push_back(ParamUI("DomainZ_________", 1.0));
    D.UI.push_back(ParamUI("RuleMinMass_____", 1.0));
    D.UI.push_back(ParamUI("RuleMaxMass_____", 1.0));
    D.UI.push_back(ParamUI("RuleMinRadius___", 0.01));
    D.UI.push_back(ParamUI("RuleMaxRadius___", 0.0));
    D.UI.push_back(ParamUI("RuleMinAmpli____", -1.0));
    D.UI.push_back(ParamUI("RuleMaxAmpli____", 1.0));
    D.UI.push_back(ParamUI("RuleMinReach____", 0.1));
    D.UI.push_back(ParamUI("RuleMaxReach____", 0.0));
    D.UI.push_back(ParamUI("RuleAdhReach____", 0.1));
    D.UI.push_back(ParamUI("ModeAdhesion____", 1));
    D.UI.push_back(ParamUI("ForceColli______", 4.0));
    D.UI.push_back(ParamUI("ForcePartic_____", 0.2));
    D.UI.push_back(ParamUI("ForceBoundary___", 0.5));
    D.UI.push_back(ParamUI("ForceDamping____", 0.5));
    D.UI.push_back(ParamUI("ForceGravity____", 0.0));
    D.UI.push_back(ParamUI("ForceAdhesion___", 0.1));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));

    D.displayModeLabel[1]= "Partic";
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool ParticLifeOrga::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool ParticLifeOrga::CheckRefresh() {
  if (D.UI[NbTypes_________].hasChanged()) isRefreshed= false;
  if (D.UI[NbPartic________].hasChanged()) isRefreshed= false;
  if (D.UI[DomainX_________].hasChanged()) isRefreshed= false;
  if (D.UI[DomainY_________].hasChanged()) isRefreshed= false;
  if (D.UI[DomainZ_________].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void ParticLifeOrga::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;
}


// Refresh the project
void ParticLifeOrga::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Get domain dimensions
  D.boxMin= {0.5 - 0.5 * D.UI[DomainX_________].D(), 0.5 - 0.5 * D.UI[DomainY_________].D(), 0.5 - 0.5 * D.UI[DomainZ_________].D()};
  D.boxMax= {0.5 + 0.5 * D.UI[DomainX_________].D(), 0.5 + 0.5 * D.UI[DomainY_________].D(), 0.5 + 0.5 * D.UI[DomainZ_________].D()};

  // Generate rules and particles if needed
  if ((int)TypeMass.size() != std::max(D.UI[NbTypes_________].I(), 1)) {
    GenerateParticleRules();
    GenerateParticleCloud();
  }

  // Generate cloud if needed
  if ((int)Type.size() != std::max(D.UI[NbPartic________].I(), 1)) {
    GenerateParticleCloud();
  }
}


// Handle UI parameter change
void ParticLifeOrga::ParamChange() {
}


// Handle keypress
void ParticLifeOrga::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  if (D.keyLetterUpperCase == 'R')
    GenerateParticleRules();
  if (D.keyLetterUpperCase == 'P')
    GenerateParticleCloud();
}


// Handle mouse action
void ParticLifeOrga::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void ParticLifeOrga::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Get and check particle system sizes
  const int nbTypes= (int)TypeMass.size();
  const int nbParticles= (int)Type.size();
  if (nbTypes <= 0) return;
  if (nbParticles <= 0) return;

  // Project to 2D
  if (D.UI[Use2DProj_______].I() > 0) {
    for (int k= 0; k < nbParticles; k++) {
      Pos[k][0]= 0.5;
      Vel[k][0]= 0.0;
    }
  }

  // Step forward simulation with explicit numerical integration
  double dt= D.UI[TimeStep________].D();
  if (D.UI[NumIntegMode____].I() == 0) {
    // Evaluate net forces acting on particles
    ComputeForces();
    // Semi-Implicit Euler integration
    for (int k= 0; k < nbParticles; k++) {
      Acc[k]= For[k] / TypeMass[Type[k]];  // at+1 = ft / m
      Vel[k]+= Acc[k] * dt;                // vt+1 = at + at+1 * dt
      Pos[k]+= Vel[k] * dt;                // xt+1 = vt + vt+1 * dt
    }
  }
  else {
    // Velocity Verlet integration - position update
    for (int k= 0; k < nbParticles; k++)
      Pos[k]+= Vel[k] * dt + 0.5 * Acc[k] * dt * dt;  // xt+1 = xt + vt * dt + 0.5 * at * dt * dt
    // Evaluate net forces acting on particles
    ComputeForces();
    // Velocity Verlet integration - acceleration and velocity update
    for (int k= 0; k < nbParticles; k++) {
      Vec::Vec3<double> oldAcc= Acc[k];
      Acc[k]= For[k] / TypeMass[Type[k]];     // at+1 = ft+1 / m
      Vel[k]+= 0.5 * (oldAcc + Acc[k]) * dt;  // vt+1 = vt + 0.5 * (at + at+1) * dt
    }
  }
}


// Draw the project
void ParticLifeOrga::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Get and check particle system sizes
  const int nbTypes= (int)TypeMass.size();
  const int nbParticles= (int)Type.size();
  if (nbTypes <= 0) return;
  if (nbParticles <= 0) return;

  // Display particles
  if (D.displayMode[1]) {
    glEnable(GL_LIGHTING);
    for (int k= 0; k < nbParticles; k++) {
      // Set particle color
      float r= 0.5, g= 0.5, b= 0.5;
      if (D.UI[ColorMode_______].I() == 0) {
        if (nbTypes > 1) Colormap::RatioToJetBrightSmooth(double(Type[k]) / (double(nbTypes) - 1.0), r, g, b);
      }
      if (D.UI[ColorMode_______].I() == 1) {
        r= 10.0 * Vel[k][0] + 0.5;
        g= 10.0 * Vel[k][1] + 0.5;
        b= 10.0 * Vel[k][2] + 0.5;
      }
      if (D.UI[ColorMode_______].I() == 2) Colormap::RatioToJetBrightSmooth(double(For[k].norm() * 10.0), r, g, b);
      glColor3f(r, g, b);
      // Draw sphere
      glPushMatrix();
      glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
      glScalef(TypeRadius[Type[k]], TypeRadius[Type[k]], TypeRadius[Type[k]]);
      glutSolidSphere(1.0, 16, 8);
      glPopMatrix();
    }
    glDisable(GL_LIGHTING);
  }
}


void ParticLifeOrga::GenerateParticleRules() {
  // Get parameters
  const int nbTypes= std::max(D.UI[NbTypes_________].I(), 1);

  // Allocate
  TypeMass= std::vector<double>(nbTypes, 0.0);
  TypeRadius= std::vector<double>(nbTypes, 0.0);
  TypeAmpli= Field::AllocNested2(nbTypes, nbTypes, 0.0);
  TypeReach= Field::AllocNested2(nbTypes, nbTypes, 0.0);

  // Generate random rule values in selected ranges
  for (int k0= 0; k0 < nbTypes; k0++) {
    TypeMass[k0]= D.UI[RuleMinMass_____].D();
    TypeRadius[k0]= D.UI[RuleMinRadius___].D();
    for (int k1= 0; k1 < nbTypes; k1++) {
      TypeAmpli[k0][k1]= D.UI[RuleMinAmpli____].D();
      TypeReach[k0][k1]= D.UI[RuleMinReach____].D();
    }

    if (D.UI[RuleMinMass_____].D() < D.UI[RuleMaxMass_____].D()) TypeMass[k0]= Random::Val(D.UI[RuleMinMass_____].D(), D.UI[RuleMaxMass_____].D());
    if (D.UI[RuleMinRadius___].D() < D.UI[RuleMaxRadius___].D()) TypeRadius[k0]= Random::Val(D.UI[RuleMinRadius___].D(), D.UI[RuleMaxRadius___].D());
    for (int k1= 0; k1 < nbTypes; k1++) {
      if (D.UI[RuleMinAmpli____].D() < D.UI[RuleMaxAmpli____].D()) TypeAmpli[k0][k1]= Random::Val(D.UI[RuleMinAmpli____].D(), D.UI[RuleMaxAmpli____].D());
      if (D.UI[RuleMinReach____].D() < D.UI[RuleMaxReach____].D()) TypeReach[k0][k1]= Random::Val(D.UI[RuleMinReach____].D(), D.UI[RuleMaxReach____].D());
    }
  }

  if (D.UI[VerboseLevel____].I() > 0) {
    printf("Mass:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      printf(" %f", TypeMass[k0]);
    }
    printf("\n");
    printf("Radius:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      printf(" %f", TypeRadius[k0]);
    }
    printf("\n");
    printf("Ampli:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      for (int k1= 0; k1 < nbTypes; k1++) {
        printf(" %f", TypeAmpli[k0][k1]);
      }
      printf("\n");
    }
    printf("Reach:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      for (int k1= 0; k1 < nbTypes; k1++) {
        printf(" %f", TypeReach[k0][k1]);
      }
      printf("\n");
    }
  }
}


void ParticLifeOrga::GenerateParticleCloud() {
  // Get parameters
  const int nbTypes= int(TypeMass.size());
  const int nbParticles= std::max(D.UI[NbPartic________].I(), 1);

  // Allocate
  Type= std::vector<int>(nbParticles);
  Pos= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));
  Vel= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));
  Acc= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));
  For= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));

  // Set random positions, velocities and types
  for (int k= 0; k < nbParticles; k++) {
    Type[k]= Random::Val(0, nbTypes - 1);
    Pos[k][0]= D.boxMin[0] + Random::Val(0.0, 1.0) * (D.boxMax[0] - D.boxMin[0]);
    Pos[k][1]= D.boxMin[1] + Random::Val(0.0, 1.0) * (D.boxMax[1] - D.boxMin[1]);
    Pos[k][2]= D.boxMin[2] + Random::Val(0.0, 1.0) * (D.boxMax[2] - D.boxMin[2]);
    // Vel[k]= {0.0, 0.0, 0.0};
    Vel[k][0]= 0.1 * (2.0 * Random::Val(0.0, 1.0) - 1.0);
    Vel[k][1]= 0.1 * (2.0 * Random::Val(0.0, 1.0) - 1.0);
    Vel[k][2]= 0.1 * (2.0 * Random::Val(0.0, 1.0) - 1.0);
    Acc[k]= {0.0, 0.0, 0.0};
    For[k]= {0.0, 0.0, 0.0};
  }

  // Project to 2D
  if (D.UI[Use2DProj_______].I() > 0) {
    for (int k= 0; k < nbParticles; k++) {
      Pos[k][0]= 0.5;
      Vel[k][0]= 0.0;
    }
  }
}


void ParticLifeOrga::ComputeForces() {
  // Get and check particle system sizes
  const int nbTypes= (int)TypeMass.size();
  const int nbParticles= (int)Type.size();
  if (nbTypes <= 0) return;
  if (nbParticles <= 0) return;

  // Reset forces
  std::fill(For.begin(), For.end(), Vec::Vec3<double>(0.0, 0.0, 0.0));

  // Compute particle forces
  for (int k0= 0; k0 < nbParticles; k0++) {
    // Interaction forces
    for (int k1= 0; k1 < nbParticles; k1++) {
      if (k0 == k1) continue;
      int typeK= Type[k0];
      int typeN= Type[k1];
      Vec::Vec3<double> distVec= Pos[k0] - Pos[k1];
      double distVal= distVec.norm();
      if (distVal == 0.0) continue;

      // Anti collision force
      double distCol= TypeRadius[typeK] + TypeRadius[typeN];
      if (distVal < distCol) {
        // double distFact= 1.0-distVal/distCol; // Linear decay
        double distFact= std::pow((distVal - distCol) / distCol, 2.0);  // Quadratic decay
        For[k0]+= D.UI[ForceColli______].D() * distFact * distVec / distVal;
      }

      // Particle attraction and repulsion force
      double distForce= TypeReach[typeN][typeK];
      // if (distVal > distCol && distVal < distForce) {
      //   // double distFact= 1.0 - std::abs((distForce + distCol - 2.0 * distVal) / (distForce - distCol));  // Linear decay
      //   double distFact= (16 * std::pow(distCol - distVal, 2.0) * std::pow(distForce - distVal, 2.0)) / std::pow(distCol - distForce, 4.0);  // Power 4 smooth decay

      if (distVal < distForce) {
        // double distFact= 1.0-distVal/distForce; // Linear decay
        double distFact= std::pow((distVal - distForce) / distForce, 2.0);  // Quadratic decay
        For[k0]-= D.UI[ForcePartic_____].D() * D.UI[ForcePartic_____].D() * TypeAmpli[typeN][typeK] * distFact * distVec / distVal;
      }

      // Adhesion force
      double distAdh= D.UI[RuleAdhReach____].D();
      if (distVal < distAdh) {
        // double distFact= 1.0; // Constant force
        // double distFact= 1.0-distVal/distAdh; // Linear decay
        double distFact= std::pow((distVal - distAdh) / distAdh, 2.0);  // Quadratic decay
        Vec::Vec3<double> diffVel= Vel[k1] - Vel[k0];
        if (D.UI[ModeAdhesion____].I() == 0) {
          diffVel.normalize();
        }
        else if (D.UI[ModeAdhesion____].I() == 1) {
          diffVel= diffVel.normalized() * std::min(Vel[k0].norm(), Vel[k1].norm());
        }
        else if (D.UI[ModeAdhesion____].I() == 2) {
          if (diffVel.norm() > std::min(Vel[k0].norm(), Vel[k1].norm())) {
            diffVel= diffVel.normalized() * std::min(Vel[k0].norm(), Vel[k1].norm());
          }
        }
        For[k0]+= D.UI[ForceAdhesion___].D() * distFact * diffVel;
      }
    }

    // Boundary repulsion forces
    double distCol= TypeRadius[Type[k0]];
    for (int dim= 0; dim < 3; dim++) {
      if (Pos[k0][dim] - distCol < D.boxMin[dim]) {
        double distFact= std::abs((Pos[k0][dim] - distCol - D.boxMin[dim]) / distCol);  // Linear repulsion
        // double distFact= std::pow((Pos[k0][dim]-distCol)/distCol, 2.0); // Quadratic repulsion
        For[k0][dim]+= D.UI[ForceBoundary___].D() * distFact;
      }
      if (Pos[k0][dim] + distCol > D.boxMax[dim]) {
        double distFact= std::abs((Pos[k0][dim] + distCol - D.boxMax[dim]) / distCol);
        // double distFact= std::pow((Pos[k0][dim]-(1.0-distCol))/distCol, 2.0);
        For[k0][dim]-= D.UI[ForceBoundary___].D() * distFact;
      }
    }

    // Gravity force toward z-
    For[k0][2]+= -0.1 * D.UI[ForceGravity____].D() * TypeMass[Type[k0]];

    // Velocity damping, i.e. friction force
    For[k0]-= D.UI[ForceDamping____].D() * Vel[k0];
  }
}
