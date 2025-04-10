#include "PosiBasedDynam.hpp"


// Standard lib
#include <cmath>
#include <vector>
#include <random>

// Algo headers
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void PosiBasedDynam::TimeIntegrate() {
  // Get parameters
  const float dt= D.UI[TimeStep________].F() / float(D.UI[NbSubSteps______].I());
  const Vec::Vec3<float> vecGrav(0.0f, 0.0f, D.UI[ForceGrav_______].F());
  const float dragValInv= 1.0f / (1.0f + dt * D.UI[ForceDrag_______].F());
  
  // Simulate for the chosen number of substeps
  for (int idxSubstep= 0; idxSubstep < D.UI[NbSubSteps______].I(); idxSubstep++) {
    // Apply external forces to velocity
    for (int k0= 0; k0 < (int)Pos.size(); k0++) {
      Vel[k0]= Vel[k0] + dt * vecGrav;      // Apply gravity
      Vel[k0]= Vel[k0] * dragValInv;        // Apply implicit drag force
    }

    // Apply random position perturbation to counteract collapse
    if (D.UI[RandPerturb_____].F() > 0.0f) {
      std::default_random_engine rng(0);
      std::uniform_real_distribution<float> distrib(-D.UI[RandPerturb_____].F(), D.UI[RandPerturb_____].F());
      for (int k0= 0; k0 < (int)Pos.size(); k0++) {
        Pos[k0][0]+= distrib(rng);
        Pos[k0][1]+= distrib(rng);
        Pos[k0][2]+= distrib(rng);
      }
    }

    // Save old positions
    PosOld= Pos;

    // Update new positions
    for (int k0= 0; k0 < (int)Pos.size(); k0++) {
      Pos[k0]= Pos[k0] + dt * Vel[k0];
    }

    // Apply position constraints
    if (D.UI[StiffBox________].F() > 0.0f) ApplyBoxDomainConstraint(D.UI[StiffBox________].F(), dt);
    if (D.UI[StiffCollision__].F() > 0.0f) ApplyNodeCollisionConstraint(D.UI[StiffCollision__].F(), dt);
    if (D.UI[StiffEdgeLength_].F() > 0.0f) ApplyEdgeLengthConstraint(D.UI[StiffEdgeLength_].F(), dt);
    if (D.UI[StiffTriArea____].F() > 0.0f) ApplyTriangleAreaConstraint(D.UI[StiffTriArea____].F(), dt);
    if (D.UI[StiffTetVolume__].F() > 0.0f) ApplyTetrahedronVolumeConstraint(D.UI[StiffTetVolume__].F(), dt);
    if (D.UI[StiffMouseBall__].F() > 0.0f) ApplyMouseBallConstraint(D.UI[StiffMouseBall__].F(), dt);

    // Deduce velocities from old and new positions
    for (int k0= 0; k0 < (int)Pos.size(); k0++) {
      Vel[k0]= (Pos[k0] - PosOld[k0]) / dt;
    }
  }
}
