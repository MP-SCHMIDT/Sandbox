#include "ParticLifeOrga.hpp"


// Standard lib
#include <cmath>
#include <vector>

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
ParticLifeOrga::ParticLifeOrga() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void ParticLifeOrga::SetActiveProject() {
  if (!isActivProj) {
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
    D.UI.push_back(ParamUI("ForceColli______", 4.0));
    D.UI.push_back(ParamUI("ForcePartic_____", 0.2));
    D.UI.push_back(ParamUI("ForceAdhesion___", 0.1));
    D.UI.push_back(ParamUI("ForceBoundary___", 0.5));
    D.UI.push_back(ParamUI("ForceDamping____", 0.5));
    D.UI.push_back(ParamUI("ForceGravity____", 0.0));
    D.UI.push_back(ParamUI("ForceAdhMode____", 1));
    D.UI.push_back(ParamUI("ForceAdhCoeff___", 0.1));
    D.UI.push_back(ParamUI("ColorMode_______", 0));
    D.UI.push_back(ParamUI("VerboseLevel____", 0.0));
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
  if ((int)ParticleTypeMass.size() != std::max(D.UI[NbTypes_________].I(), 1)) {
    GenerateParticleRules();
    GenerateParticleCloud();
  }

  // Generate cloud if needed
  if ((int)ParticleType.size() != std::max(D.UI[NbPartic________].I(), 1)) {
    GenerateParticleCloud();
  }
}


// Handle keypress
void ParticLifeOrga::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable

  if (key == 'R')
    GenerateParticleRules();
  if (key == 'P')
    GenerateParticleCloud();
}


// Animate the project
void ParticLifeOrga::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Get and check particle system sizes
  const int nbTypes= (int)ParticleTypeMass.size();
  const int nbParticles= (int)ParticleType.size();
  if (nbTypes <= 0) return;
  if (nbParticles <= 0) return;

  // Project to 2D
  if (D.UI[Use2DProj_______].B()) {
    for (int k= 0; k < nbParticles; k++) {
      ParticlePos[k][0]= 0.5;
      ParticleVel[k][0]= 0.0;
    }
  }

  // Step forward simulation with explicit numerical integration
  if (D.UI[NumIntegMode____].I() == 0) {
    // Evaluate net forces acting on particles
    ComputeForces();
    // Euler integration
    for (int k= 0; k < nbParticles; k++) {
      ParticleAcc[k]= ParticleFor[k] / ParticleTypeMass[ParticleType[k]];  // at+1 = ft / m
      ParticleVel[k]+= ParticleAcc[k] * D.UI[TimeStep________].D();        // vt+1 = at + at+1 * dt
      ParticlePos[k]+= ParticleVel[k] * D.UI[TimeStep________].D();        // xt+1 = vt + vt+1 * dt
    }
  }
  else {
    // Velocity Verlet integration - position update
    for (int k= 0; k < nbParticles; k++)
      ParticlePos[k]+= ParticleVel[k] * D.UI[TimeStep________].D() + 0.5 * ParticleAcc[k] * D.UI[TimeStep________].D() * D.UI[TimeStep________].D();  // xt+1 = xt + vt * dt + 0.5 * at * dt * dt
    // Evaluate net forces acting on particles
    ComputeForces();
    // Velocity Verlet integration - acceleration and velocity update
    for (int k= 0; k < nbParticles; k++) {
      Vec::Vec3<double> oldAcc= ParticleAcc[k];
      ParticleAcc[k]= ParticleFor[k] / ParticleTypeMass[ParticleType[k]];             // at+1 = ft+1 / m
      ParticleVel[k]+= 0.5 * (oldAcc + ParticleAcc[k]) * D.UI[TimeStep________].D();  // vt+1 = vt + 0.5 * (at + at+1) * dt
    }
  }

  // // Update density trail
  // // Get field sizes and init if non existant
  // int nbX, nbY, nbZ;
  // SrtUtil::GetFieldDimensions(densityField, nbX, nbY, nbZ);
  // if (int(displacementField.size()) != nbX || int(displacementField[0].size()) != nbY || int(displacementField[0][0].size()) != nbZ) {
  //   displacementField= vector<vector<vector<Vec::Vec3<double>>>>(nbX, vector<vector<Vec::Vec3<double>>>(nbY, vector<Vec::Vec3<double>>(nbZ, {0.0, 0.0, 0.0})));
  // }
  // double stepX, stepY, stepZ, voxDiag, startX, startY, startZ;
  // SrtUtil::GetVoxelSizes(nbX, nbY, nbZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ, voxDiag);
  // SrtUtil::GetVoxelStart(D.boxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  // // Loop through particles and leave trail
  // for (int k= 0; k < nbParticles; k++) {
  //   int x= std::min(std::max(int(std::round((ParticlePos[k][0] - D.boxMin[0]) / stepX)), 0), nbX - 1);
  //   int y= std::min(std::max(int(std::round((ParticlePos[k][1] - D.boxMin[1]) / stepY)), 0), nbY - 1);
  //   int z= std::min(std::max(int(std::round((ParticlePos[k][2] - D.boxMin[2]) / stepZ)), 0), nbZ - 1);
  //   densityField[x][y][z]= std::min(densityField[x][y][z] + D.UI[plconRatio5_____].D() * D.UI[plconRatio5_____].D(), 1.0);
  //   displacementField[x][y][z]= displacementField[x][y][z] + D.UI[plconRatio5_____].D() * D.UI[plconRatio5_____].D() * ParticleVel[k];
  // }
  // SrtAlgoField::SmoothScalarField(densityField);
  // SrtAlgoField::SmoothVectorField(displacementField);
  // // Trail map decay
  // for (int x= 0; x < nbX; x++) {
  //   for (int y= 0; y < nbY; y++) {
  //     for (int z= 0; z < nbZ; z++) {
  //       densityField[x][y][z]= std::max(densityField[x][y][z] - D.UI[plconRatio4_____].D() * D.UI[plconRatio4_____].D(), 0.0);
  //       displacementField[x][y][z]= displacementField[x][y][z] * (1.0 - D.UI[plconRatio4_____].D() * D.UI[plconRatio4_____].D());
  //     }
  //   }
  // }

  // // Update density trail isosurface
  // bool showIso= true;
  // if (showIso) {
  //   static vector<vector<vector<double>>> tmpDensityField;
  //   if (int(tmpDensityField.size()) != nbX + 2 || int(tmpDensityField[0].size()) != nbY + 2 || int(tmpDensityField[0][0].size()) != nbZ + 2)
  //     tmpDensityField= vector<vector<vector<double>>>(nbX + 2, vector<vector<double>>(nbY + 2, vector<double>(nbZ + 2, 0.0)));
  //   array<double, 3> expandedBBoxMin;
  //   array<double, 3> expandedBBoxMax;
  //   double scaleX= (D.boxMax[0] - D.boxMin[0]) / double(nbX);
  //   double scaleY= (D.boxMax[1] - D.boxMin[1]) / double(nbY);
  //   double scaleZ= (D.boxMax[2] - D.boxMin[2]) / double(nbZ);
  //   expandedBBoxMin= {D.boxMin[0] - scaleX, D.boxMin[1] - scaleY, D.boxMin[2] - scaleZ};
  //   expandedBBoxMax= {D.boxMax[0] + scaleX, D.boxMax[1] + scaleY, D.boxMax[2] + scaleZ};
  //   for (int x= 0; x < nbX; x++)
  //     for (int y= 0; y < nbY; y++)
  //       for (int z= 0; z < nbZ; z++)
  //         tmpDensityField[x + 1][y + 1][z + 1]= densityField[x][y][z];
  //   SrtMarchingCubes::ComputeMarchingCubes(tmpDensityField, expandedBBoxMin, expandedBBoxMax, D.UI[DensiCut________].D(), vertsIso, trisIso);
  //   normalsIso.clear();
  // }
}


// Draw the project
void ParticLifeOrga::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Get and check particle system sizes
  const int nbTypes= (int)ParticleTypeMass.size();
  const int nbParticles= (int)ParticleType.size();
  if (nbTypes <= 0) return;
  if (nbParticles <= 0) return;

  // Display particles
  if (D.displayMode1) {
    glEnable(GL_LIGHTING);
    for (int k= 0; k < nbParticles; k++) {
      // Set particle color
      float r= 0.5, g= 0.5, b= 0.5;
      if (D.UI[ColorMode_______].I() == 0) {
        if (nbTypes > 1) Colormap::RatioToJetBrightSmooth(double(ParticleType[k]) / (double(nbTypes) - 1.0), r, g, b);
      }
      if (D.UI[ColorMode_______].I() == 1) {
        r= 10.0 * ParticleVel[k][0] + 0.5;
        g= 10.0 * ParticleVel[k][1] + 0.5;
        b= 10.0 * ParticleVel[k][2] + 0.5;
      }
      if (D.UI[ColorMode_______].I() == 2) Colormap::RatioToJetBrightSmooth(double(ParticleFor[k].norm() * 10.0), r, g, b);
      glColor3f(r, g, b);
      // Draw sphere
      glPushMatrix();
      glTranslatef(ParticlePos[k][0], ParticlePos[k][1], ParticlePos[k][2]);
      glScalef(ParticleTypeRadius[ParticleType[k]], ParticleTypeRadius[ParticleType[k]], ParticleTypeRadius[ParticleType[k]]);
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
  ParticleTypeMass= std::vector<double>(nbTypes, 0.0);
  ParticleTypeRadius= std::vector<double>(nbTypes, 0.0);
  ParticleTypeAmpli= Field::AllocField2D(nbTypes, nbTypes, 0.0);
  ParticleTypeReach= Field::AllocField2D(nbTypes, nbTypes, 0.0);

  // Generate random rule values in selected ranges
  for (int k0= 0; k0 < nbTypes; k0++) {
    ParticleTypeMass[k0]= D.UI[RuleMinMass_____].D();
    ParticleTypeRadius[k0]= D.UI[RuleMinRadius___].D();
    for (int k1= 0; k1 < nbTypes; k1++) {
      ParticleTypeAmpli[k0][k1]= D.UI[RuleMinAmpli____].D();
      ParticleTypeReach[k0][k1]= D.UI[RuleMinReach____].D();
    }

    if (D.UI[RuleMinMass_____].D() < D.UI[RuleMaxMass_____].D()) ParticleTypeMass[k0]= Random::Val(D.UI[RuleMinMass_____].D(), D.UI[RuleMaxMass_____].D());
    if (D.UI[RuleMinRadius___].D() < D.UI[RuleMaxRadius___].D()) ParticleTypeRadius[k0]= Random::Val(D.UI[RuleMinRadius___].D(), D.UI[RuleMaxRadius___].D());
    for (int k1= 0; k1 < nbTypes; k1++) {
      if (D.UI[RuleMinAmpli____].D() < D.UI[RuleMaxAmpli____].D()) ParticleTypeAmpli[k0][k1]= Random::Val(D.UI[RuleMinAmpli____].D(), D.UI[RuleMaxAmpli____].D());
      if (D.UI[RuleMinReach____].D() < D.UI[RuleMaxReach____].D()) ParticleTypeReach[k0][k1]= Random::Val(D.UI[RuleMinReach____].D(), D.UI[RuleMaxReach____].D());
    }
  }

  if (D.UI[VerboseLevel____].B()) {
    printf("Mass:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      printf(" %f", ParticleTypeMass[k0]);
    }
    printf("\n");
    printf("Radius:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      printf(" %f", ParticleTypeRadius[k0]);
    }
    printf("\n");
    printf("Ampli:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      for (int k1= 0; k1 < nbTypes; k1++) {
        printf(" %f", ParticleTypeAmpli[k0][k1]);
      }
      printf("\n");
    }
    printf("Reach:\n");
    for (int k0= 0; k0 < nbTypes; k0++) {
      for (int k1= 0; k1 < nbTypes; k1++) {
        printf(" %f", ParticleTypeReach[k0][k1]);
      }
      printf("\n");
    }
  }
}


void ParticLifeOrga::GenerateParticleCloud() {
  // Get parameters
  const int nbTypes= int(ParticleTypeMass.size());
  const int nbParticles= std::max(D.UI[NbPartic________].I(), 1);

  // Allocate
  ParticleType= std::vector<int>(nbParticles);
  ParticlePos= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));
  ParticleVel= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));
  ParticleAcc= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));
  ParticleFor= std::vector<Vec::Vec3<double>>(nbParticles, Vec::Vec3<double>(0.0, 0.0, 0.0));

  // Set random positions, velocities and types
  for (int k= 0; k < nbParticles; k++) {
    ParticleType[k]= Random::Val(0, nbTypes - 1);
    ParticlePos[k][0]= D.boxMin[0] + Random::Val(0.0, 1.0) * (D.boxMax[0] - D.boxMin[0]);
    ParticlePos[k][1]= D.boxMin[1] + Random::Val(0.0, 1.0) * (D.boxMax[1] - D.boxMin[1]);
    ParticlePos[k][2]= D.boxMin[2] + Random::Val(0.0, 1.0) * (D.boxMax[2] - D.boxMin[2]);
    ParticleVel[k][0]= 0.1 * (2.0 * Random::Val(0.0, 1.0) - 1.0);
    ParticleVel[k][1]= 0.1 * (2.0 * Random::Val(0.0, 1.0) - 1.0);
    ParticleVel[k][2]= 0.1 * (2.0 * Random::Val(0.0, 1.0) - 1.0);
    ParticleAcc[k]= {0.0, 0.0, 0.0};
    ParticleFor[k]= {0.0, 0.0, 0.0};
  }

  // Project to 2D
  if (D.UI[Use2DProj_______].B()) {
    for (int k= 0; k < nbParticles; k++) {
      ParticlePos[k][0]= 0.5;
      ParticleVel[k][0]= 0.0;
    }
  }
}


void ParticLifeOrga::ComputeForces() {
  // Get and check particle system sizes
  const int nbTypes= (int)ParticleTypeMass.size();
  const int nbParticles= (int)ParticleType.size();
  if (nbTypes <= 0) return;
  if (nbParticles <= 0) return;

  // Reset forces and neighbors counts
  std::fill(ParticleFor.begin(), ParticleFor.end(), Vec::Vec3<double>(0.0, 0.0, 0.0));

  // Compute particle forces
  for (int k0= 0; k0 < nbParticles; k0++) {
    // Interaction forces
    for (int k1= 0; k1 < nbParticles; k1++) {
      if (k0 == k1) continue;
      int typeK= ParticleType[k0];
      int typeN= ParticleType[k1];
      Vec::Vec3<double> vec1to0= ParticlePos[k0] - ParticlePos[k1];
      double dist= vec1to0.norm();
      vec1to0= vec1to0 / dist;

      // Anti collision force
      double distCol= ParticleTypeRadius[typeK] + ParticleTypeRadius[typeN];
      if (dist < distCol) {
        // double distFact= 1.0-dist/distCol; // Linear decay
        double distFact= std::pow((dist - distCol) / distCol, 2.0);  // Quadratic decay
        ParticleFor[k0]+= D.UI[ForceColli______].D() * distFact * vec1to0;
      }

      // Particle attraction and repulsion force
      double distForce= ParticleTypeReach[typeN][typeK];
      if (dist > distCol && dist < distForce) {
        // double distFact= 1.0-std::abs((distForce+distCol-2.0*dist)/(distForce-distCol)); // Linear decay
        double distFact= (16 * std::pow(distCol - dist, 2.0) * std::pow(distForce - dist, 2.0)) / std::pow(distCol - distForce, 4.0);  // Power 4 smooth decay
        ParticleFor[k0]-= D.UI[ForcePartic_____].D() * D.UI[ForcePartic_____].D() * ParticleTypeAmpli[typeN][typeK] * distFact * vec1to0;
      }

      // Adhesion force
      double distAdh= D.UI[ForceAdhesion___].D();
      if (dist < distAdh) {
        // double distFact= 1.0; // Constant force
        // double distFact= 1.0-dist/distAdh; // Linear decay
        double distFact= std::pow((dist - distAdh) / distAdh, 2.0);  // Quadratic decay
        Vec::Vec3<double> diffVel= ParticleVel[k1] - ParticleVel[k0];
        if (D.UI[ForceAdhMode____].I() == 0) {
          diffVel.normalize();
        }
        else if (D.UI[ForceAdhMode____].I() == 1) {
          diffVel= diffVel.normalized() * std::min(ParticleVel[k0].norm(), ParticleVel[k1].norm());
        }
        else if (D.UI[ForceAdhMode____].I() == 2) {
          if (diffVel.norm() > std::min(ParticleVel[k0].norm(), ParticleVel[k1].norm())) {
            diffVel= diffVel.normalized() * std::min(ParticleVel[k0].norm(), ParticleVel[k1].norm());
          }
        }
        ParticleFor[k0]+= D.UI[ForceAdhCoeff___].D() * distFact * diffVel;
      }
    }

    // Boundary repulsion forces
    double distCol= ParticleTypeRadius[ParticleType[k0]];
    for (int dim= 0; dim < 3; dim++) {
      if (ParticlePos[k0][dim] - distCol < D.boxMin[dim]) {
        double distFact= std::abs((ParticlePos[k0][dim] - distCol - D.boxMin[dim]) / distCol);
        // double distFact= std::pow((ParticlePos[k0][dim]-distCol)/distCol, 2.0);
        ParticleFor[k0][dim]+= D.UI[ForceBoundary___].D() * distFact;
      }
      if (ParticlePos[k0][dim] + distCol > D.boxMax[dim]) {
        double distFact= std::abs((ParticlePos[k0][dim] + distCol - D.boxMax[dim]) / distCol);
        // double distFact= std::pow((ParticlePos[k0][dim]-(1.0-distCol))/distCol, 2.0);
        ParticleFor[k0][dim]-= D.UI[ForceBoundary___].D() * distFact;
      }
    }

    // Gravity force toward z-
    ParticleFor[k0][2]+= -0.1 * D.UI[ForceGravity____].D() * ParticleTypeMass[ParticleType[k0]];

    // Velocity damping, i.e. friction force
    ParticleFor[k0]-= D.UI[ForceDamping____].D() * ParticleVel[k0];
  }
}
