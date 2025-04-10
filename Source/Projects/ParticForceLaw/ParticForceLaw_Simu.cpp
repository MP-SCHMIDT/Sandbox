#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <numbers>
#include <vector>

// Algo headers
#include "Geom/BoxGrid.hpp"
#include "Geom/MarchingCubes.hpp"
#include "Geom/SpatialBucket3D.hpp"
#include "Type/Field.hpp"
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void ParticForceLaw::AddForceInteractions(const std::vector<Vec::Vec3<float>>& iPos,
                                          const std::vector<Vec::Vec3<float>>& iVel,
                                          std::vector<Vec::Vec3<float>>& oFor) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset vectors
  std::fill(ForceMag.begin(), ForceMag.end(), 0.0f);
  std::fill(Neighbors.begin(), Neighbors.end(), 0);

  // Precompute values
  const float surfArea= 4.0f * std::numbers::pi * D.UI[LatticePitch____].F() * D.UI[LatticePitch____].F();
  float maxForceLawRange= 0.0f;
  std::vector<float> ForceLawRangesSqr= ForceLawRange;
  for (int idxMat= 0; idxMat < (int)ForceLaw.size(); idxMat++) {
    maxForceLawRange= std::max(maxForceLawRange, ForceLawRange[idxMat]);
    ForceLawRangesSqr[idxMat]= ForceLawRange[idxMat] * ForceLawRange[idxMat];
  }

  // Allocate and fill the spatial hash
  Vec::Vec3<float> boxMin= {(float)D.boxMin[0], (float)D.boxMin[1], (float)D.boxMin[2]};
  Vec::Vec3<float> boxMax= {(float)D.boxMax[0], (float)D.boxMax[1], (float)D.boxMax[2]};
  if (D.UI[BucketFitDomain_].I() > 0) { 
    boxMin= iPos[0];
    boxMax= iPos[0];
    for (Vec::Vec3<float> posPoint : iPos) {
      boxMin= boxMin.cwiseMin(posPoint);
      boxMax= boxMax.cwiseMax(posPoint);
    }
    if (D.UI[BucketFitDomain_].I() == 1) {
      boxMin= boxMin.cwiseMax({(float)D.boxMin[0], (float)D.boxMin[1], (float)D.boxMin[2]});
      boxMax= boxMax.cwiseMin({(float)D.boxMax[0], (float)D.boxMax[1], (float)D.boxMax[2]});
    }
  }
  if (!((boxMax - boxMin).minCoeff() >= 0.0f)) return;
  Bucket3D.SetDim(boxMin, boxMax, maxForceLawRange, D.UI[BucketCapacity__].I(), D.UI[BucketMaxCount__].I());
  BucketOverflown= Bucket3D.Fill(iPos);

  // Sweep through the points
  const Vec::Vec3<float> vecOffset(maxForceLawRange);
  #pragma omp parallel for
  for (int k0= 0; k0 < (int)iPos.size(); k0++) {
    // Sweep through the candidate cells withing reach of the point
    int minX, minY, minZ, maxX, maxY, maxZ;
    Bucket3D.GetCell(iPos[k0] - vecOffset, minX, minY, minZ);
    Bucket3D.GetCell(iPos[k0] + vecOffset, maxX, maxY, maxZ);
    for (int cellX= std::max(0, minX); cellX <= std::min(maxX, Bucket3D.nX - 1); cellX++) {
      for (int cellY= std::max(0, minY); cellY <= std::min(maxY, Bucket3D.nY - 1); cellY++) {
        for (int cellZ= std::max(0, minZ); cellZ <= std::min(maxZ, Bucket3D.nZ - 1); cellZ++) {
          // Test all the potential neighbor points in the candidate cells
          for (int k1 : Bucket3D.buckets[Bucket3D.GetFlatIdx(cellX, cellY, cellZ)]) {
            if (k0 == k1) continue;
            // Precompute distances
            const Vec::Vec3<float> distVec= iPos[k0] - iPos[k1];
            const float distSquared= distVec.normSquared();
            if (distSquared > ForceLawRangesSqr[Mat[k0]] && distSquared > ForceLawRangesSqr[Mat[k1]]) continue;
            const float distVal= std::sqrt(distSquared);
            const Vec::Vec3<float> distVecUnit= distVec / distVal;
            // Get linear interpolation of force laws for the given distance
            float forceVal= 0.0f;
            for (int k : {k0, k1}) {
              const float valFloat= distVal / (ForceLawStep[Mat[k]] * D.UI[LatticePitch____].F());
              const int low= std::min(std::max((int)std::floor(valFloat), 0), (int)ForceLaw[Mat[k]].size() - 1);
              const int upp= std::min(std::max(low + 1,                   0), (int)ForceLaw[Mat[k]].size() - 1);
              const float ratio= valFloat - (float)low;
              forceVal+= 0.5f * ((1.0f - ratio) * ForceLaw[Mat[k]][low] + (ratio)*ForceLaw[Mat[k]][upp]);
            }
            // Apply inter-particle force
            ForceMag[k0]+= std::abs(forceVal) * surfArea;
            oFor[k0]+= forceVal * surfArea * distVecUnit;
            // Apply inter-particle damping proportional to radial velocity
            const float radialVel= (iVel[k0] - iVel[k1]).dot(distVecUnit);
            ForceMag[k0]+= (1.0f - distVal / maxForceLawRange) * D.UI[DampingRadRel___].F() * ForceLaw[Mat[k0]][0] * std::abs(radialVel) * surfArea;
            oFor[k0]-= (1.0f - distVal / maxForceLawRange) * D.UI[DampingRadRel___].F() * ForceLaw[Mat[k0]][0] * radialVel * surfArea * distVecUnit;
            // Increment neighbor count for display
            Neighbors[k0]++;
          }
        }
      }
    }
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("ForcesT %f\n", Timer::PopTimer());
}


void ParticForceLaw::StepSimulation() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Precompute values
  const float dt= D.UI[TimeStep________].F();
  const float invMass= 1.0f / (D.UI[MaterialDensity_].F() * std::pow(D.UI[LatticePitch____].F(), 3.0f));
  const float damping= std::min(std::max(D.UI[DampingVelRel___].F(), 0.0f), 1.0f);
  const int ForceControl= D.UI[ForceControl____].I();
  const int use2D= D.UI[ConstrainDim2D__].I();
  const Vec::Vec3<float> ForNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> ForPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());
  const Vec::Vec3<float> VelNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> VelPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());
  const float surfArea= 4.0f * std::numbers::pi * D.UI[LatticePitch____].F() * D.UI[LatticePitch____].F();

  // Semi-implicit Euler integration

  // Compute forces
  std::fill(For.begin(), For.end(), Vec::Vec3<float>(0.0f, 0.0f, 0.0f)); // Reset forces
  AddForceInteractions(Pos, Vel, For);                                   // f(t)
  // Apply BC forces
  for (int k= 0; k < (int)Pos.size(); k++) {
    // Force controllers math
    // Equations of motion    xt+1 = xt + dt*vt + dt*dt*ft/m
    // Equations of motion    vt+1 = vt + dt*ft/m
    // Position controller    xt+1 == xbc  =>  ft == m*(xbc - xt - dt*vt) / dt*dt
    // Velocity controller    vt+1 == vbc  =>  ft == m*(vbc - vt) / dt
    if ((ForceControl == 1 && BCPos[k] != 0) ||
        (ForceControl == 2 && BCPos[k] != 0) ||
        (ForceControl == 2 && BCVel[k] != 0)) {  // Force controller for position and velocity using target reference position
      Vec::Vec3<float> ErrVec= Ref[k] - Pos[k];
      if (ErrVec.normSquared() > 0.0f)
        For[k]+= D.UI[BCPosCoeff______].F() * surfArea * (ErrVec.normalized() * ErrVec.normSquared() - dt * Vel[k]) / (dt * dt);
    }
    else if (ForceControl == 1 && BCVel[k] != 0) {  // Force controller for target velocity
      Vec::Vec3<float> ErrVec= ((BCVel[k] < 0) ? (VelNega) : (VelPosi)) - Vel[k];
      if (ErrVec.normSquared() > 0.0f)
        For[k]+= D.UI[BCVelCoeff______].F() * surfArea * ErrVec.normalized() * ErrVec.normSquared() / dt;
    }
    else if (BCFor[k] != 0) {  // Force controller for force
      For[k]+= (BCFor[k] < 0) ? (ForNega) : (ForPosi);
    }
  }

  // Update velocities
  for (int k= 0; k < (int)Pos.size(); k++) {
    Vel[k]+= dt * For[k] * invMass;  // v(t+1) = v(t) + Δt * f(t) / m
    Vel[k]*= (1.0f - damping);       // v(t+1)= ζ v(t+1)
    if (ForceControl == 0) {
      if (BCVel[k] != 0) Vel[k]= (BCVel[k] < 0) ? VelNega : VelPosi;  // Override velocity
      if (BCPos[k] != 0) Vel[k].set(0.0f, 0.0f, 0.0f);                // Override velocity
    }
    if (use2D > 0 && use2D <= 3) Vel[k][use2D - 1]= 0.0f;                   // Constrain velocity to 2D
    if (BCVel[k] != 0) Ref[k]+= dt * ((BCVel[k] < 0) ? VelNega : VelPosi);  // Update reference position
  }

  // Update positions
  for (int k= 0; k < (int)Pos.size(); k++) {
    Pos[k]+= dt * Vel[k];                     // x(t+1) = x(t) + Δt * v(t+1)
    if (ForceControl == 0) {
      if (BCPos[k] != 0) Pos[k]= Ref[k];     // Overwrite position
    }
  }

  // Advance time
  SimTime+= dt;

  // Scatter plot of sensor data
  if ((int)D.Scatter.size() < RunID) D.Scatter.resize(RunID);
  D.Scatter[RunID - 1].name= "ForceDisp";
  float sumPos= 0.0f;
  float sumFor= 0.0f;
  int count= 0;
  for (int k= 0; k < (int)Pos.size(); k++) {
    if (Sensor[k]) {
      count++;
      sumPos+= Pos[k][2];
      sumFor+= For[k].norm();
    }
  }
  D.Scatter[RunID - 1].val.push_back(std::array<double, 2>{sumPos / (float)count, sumFor / (float)count});

  if (D.UI[VerboseLevel____].I() >= 1) printf("StepT %f\n", Timer::PopTimer());
}


void ParticForceLaw::ComputeMetaballs() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  MetaballIsUpdated= true;
  Verts.clear();
  Tris.clear();

  if (D.UI[MetaballVoxSize_].F() <= 0.0f) return;

  const float metaballSize= D.UI[LatticePitch____].F();
  const int tmpNX= std::max((int)std::round((D.boxMax[0] - D.boxMin[0]) / D.UI[MetaballVoxSize_].D()), 1);
  const int tmpNY= std::max((int)std::round((D.boxMax[1] - D.boxMin[1]) / D.UI[MetaballVoxSize_].D()), 1);
  const int tmpNZ= std::max((int)std::round((D.boxMax[2] - D.boxMin[2]) / D.UI[MetaballVoxSize_].D()), 1);
  Field::Field3<double> tmpField(tmpNX, tmpNY, tmpNZ, 0.0);
  double stepX, stepY, stepZ;
  BoxGrid::GetVoxelSizes(tmpNX, tmpNY, tmpNZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ);
  for (int k= 0; k < (int)Pos.size(); k++) {
    const int idxXBeg= (int)std::floor((float)tmpNX * (Pos[k][0] - 2.0 * metaballSize - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxYBeg= (int)std::floor((float)tmpNY * (Pos[k][1] - 2.0 * metaballSize - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZBeg= (int)std::floor((float)tmpNZ * (Pos[k][2] - 2.0 * metaballSize - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    const int idxXEnd= (int)std::floor((float)tmpNX * (Pos[k][0] + 2.0 * metaballSize - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]));
    const int idxYEnd= (int)std::floor((float)tmpNY * (Pos[k][1] + 2.0 * metaballSize - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]));
    const int idxZEnd= (int)std::floor((float)tmpNZ * (Pos[k][2] + 2.0 * metaballSize - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]));
    for (int x= std::max(0, idxXBeg); x <= std::min(idxXEnd, tmpNX - 1); x++) {
      for (int y= std::max(0, idxYBeg); y <= std::min(idxYEnd, tmpNY - 1); y++) {
        for (int z= std::max(0, idxZBeg); z <= std::min(idxZEnd, tmpNZ - 1); z++) {
          const Vec::Vec3<float> pos(D.boxMin[0] + (x + 0.5f) * stepX, D.boxMin[1] + (y + 0.5f) * stepY, D.boxMin[2] + (z + 0.5f) * stepZ);
          if ((pos - Pos[k]).normSquared() < (2.0 * metaballSize) * (2.0 * metaballSize)) {
            // TODO tweak decay function for cleaner results
            tmpField.at(x, y, z)+= std::max(0.0, 1.0 - (pos - Pos[k]).norm() / metaballSize);
          }
        }
      }
    }
  }

  MarchingCubes::BuildMesh(tmpField.nX, tmpField.nY, tmpField.nZ, true, true, D.UI[MetaballIsoval__].F(), D.boxMin, D.boxMax, tmpField.data, Verts, Tris);

  if (D.UI[VerboseLevel____].I() >= 1) printf("MetaballT %f\n", Timer::PopTimer());
}
