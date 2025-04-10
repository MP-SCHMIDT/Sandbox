#include "PosiBasedDynam.hpp"


// Standard lib
#include <cmath>
#include <vector>

// Algo headers
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void PosiBasedDynam::ApplyBoxDomainConstraint(const float iStiffness, const float iTimestep) {
  const float softnessCoeff= 1.0f / (iStiffness * iTimestep * iTimestep);
  const Vec::Vec3<float> boxMin((float)D.boxMin[0], (float)D.boxMin[1], (float)D.boxMin[2]);
  const Vec::Vec3<float> boxMax((float)D.boxMax[0], (float)D.boxMax[1], (float)D.boxMax[2]);
  for (int k0= 0; k0 < (int)Pos.size(); k0++) {
    const Vec::Vec3<float> newPos= Pos[k0].cwiseMin(boxMax).cwiseMax(boxMin);
    if (newPos != Pos[k0]) {
      const float lambda= 1.0f / (MassInv[k0] + softnessCoeff);
      Pos[k0]+= lambda * MassInv[k0] * (newPos - Pos[k0]);
    }
  }
}


void PosiBasedDynam::ApplyNodeCollisionConstraint(const float iStiffness, const float iTimestep) {
  const float softnessCoeff= 1.0f / (iStiffness * iTimestep * iTimestep);
  const float targetDistance= 2.0f * D.UI[RadParticl______].F();
  const float targetDistanceSqr= targetDistance * targetDistance;
  // Nested loop to find pairwise collisions
  // O(n^2) implementation because Gauss Seidel style relaxation
  // Possible to go O(n) with spatial hash but no longer stable sequential GS
  for (int k0= 0; k0 < D.UI[NumParticl______].I(); k0++) {
    for (int k1= k0 + 1; k1 < (int)Pos.size(); k1++) {
      const float distSqr= (Pos[k0] - Pos[k1]).normSquared();
      if (distSqr > 0.0f && distSqr < targetDistanceSqr) {
        // Apply the pairwise collision position correction
        const float dist= std::sqrt(distSqr);
        const float lambda= -(dist - targetDistance) / (MassInv[k0] + MassInv[k1] + softnessCoeff);
        const Vec::Vec3<float> grad0= (Pos[k0] - Pos[k1]) / dist;
        Pos[k0]+= lambda * MassInv[k0] * grad0;
        Pos[k1]-= lambda * MassInv[k1] * grad0;
      }
    }
  }
}


void PosiBasedDynam::ApplyEdgeLengthConstraint(const float iStiffness, const float iTimestep) {
  const float softnessCoeff= 1.0f / (iStiffness * iTimestep * iTimestep);
  for (int idxEdge= 0; idxEdge < (int)Edge.size(); idxEdge++) {
    const int k0= Edge[idxEdge][0], k1= Edge[idxEdge][1];
    const float distSqr= (Pos[k0] - Pos[k1]).normSquared();
    if (distSqr > 0.0f) {
      // Apply the edge length correction
      const float dist= std::sqrt(distSqr);
      const float lambda= -(dist - EdgeLen[idxEdge]) / (MassInv[k0] + MassInv[k1] + softnessCoeff);
      const Vec::Vec3<float> grad0= (Pos[k0] - Pos[k1]) / dist;
      Pos[k0]+= lambda * MassInv[k0] * grad0;
      Pos[k1]-= lambda * MassInv[k1] * grad0;
    }
  }
}


void PosiBasedDynam::ApplyTriangleAreaConstraint(const float iStiffness, const float iTimestep) {
  const float softnessCoeff= 1.0f / (iStiffness * iTimestep * iTimestep);
  for (int idxTri= 0; idxTri < (int)Tri.size(); idxTri++) {
    const int k0= Tri[idxTri][0], k1= Tri[idxTri][1], k2= Tri[idxTri][2];
    const Vec::Vec3<float> p0= Pos[k0], p1= Pos[k1], p2= Pos[k2];
    const Vec::Vec3<float> n= ((p1-p0).cross(p2-p0)).normalized();
    const float area= GetTriArea(p0, p1, p2);
    const Vec::Vec3<float> grad0= 0.5f * n.cross(p1-p2);
    const Vec::Vec3<float> grad1= 0.5f * n.cross(p2-p0);
    const Vec::Vec3<float> grad2= 0.5f * n.cross(p0-p1);
    const float lambda= (area - TriArea[idxTri]) / (MassInv[k0]*grad0.normSquared() +
                                                    MassInv[k1]*grad1.normSquared() +
                                                    MassInv[k2]*grad2.normSquared() + softnessCoeff);
    Pos[k0]+= lambda * MassInv[k0] * grad0;
    Pos[k1]+= lambda * MassInv[k1] * grad1;
    Pos[k2]+= lambda * MassInv[k2] * grad2;
  }
}


void PosiBasedDynam::ApplyTetrahedronVolumeConstraint(const float iStiffness, const float iTimestep) {
  const float softnessCoeff= 1.0f / (iStiffness * iTimestep * iTimestep);
  for (int idxTet= 0; idxTet < (int)Tet.size(); idxTet++) {
    // Apply the tet volume correction
    const int k0= Tet[idxTet][0], k1= Tet[idxTet][1], k2= Tet[idxTet][2], k3= Tet[idxTet][3];
    const Vec::Vec3<float> p0= Pos[k0], p1= Pos[k1], p2= Pos[k2], p3= Pos[k3];
    const float vol= GetTetVolume(p0, p1, p2, p3);
    const Vec::Vec3<float> grad0= (p3-p1).cross(p2-p1);
    const Vec::Vec3<float> grad1= (p2-p0).cross(p3-p0);
    const Vec::Vec3<float> grad2= (p3-p0).cross(p1-p0);
    const Vec::Vec3<float> grad3= (p1-p0).cross(p2-p0);
    const float lambda= -6.0f * (vol - TetVol[idxTet]) / (MassInv[k0]*grad0.normSquared() +
                                                          MassInv[k1]*grad1.normSquared() +
                                                          MassInv[k2]*grad2.normSquared() +
                                                          MassInv[k3]*grad3.normSquared() + softnessCoeff);
    Pos[k0]+= lambda * MassInv[k0] * grad0;
    Pos[k1]+= lambda * MassInv[k1] * grad1;
    Pos[k2]+= lambda * MassInv[k2] * grad2;
    Pos[k3]+= lambda * MassInv[k3] * grad3;
  }
}


void PosiBasedDynam::ApplyMouseBallConstraint(const float iStiffness, const float iTimestep) {
  if (D.mouseMiddleButtonState > 0) {
    const float softnessCoeff= 1.0f / (iStiffness * iTimestep * iTimestep);
    const float targetDist= D.UI[RadMouse________].F() + D.UI[RadParticl______].F();
    const float targetDistSqr= targetDist * targetDist;
    const Vec::Vec3<float> mouse(D.mouseProjX[0], D.mouseProjX[1], D.mouseProjX[2]);
    for (int k0= 0; k0 < (int)Pos.size(); k0++) {
      const Vec::Vec3<float> vec= Pos[k0]-mouse;
      const float distSqr= vec.normSquared();
      if (distSqr > 0.0f && distSqr < targetDistSqr) {
        // Apply the repulsion position correction
        const Vec::Vec3<float> newPos= mouse + vec.normalized() * targetDist;
        const float lambda= 1.0f / (MassInv[k0] + softnessCoeff);
        Pos[k0]+= lambda * MassInv[k0] * (newPos - Pos[k0]);
      }
    }
  }
}
