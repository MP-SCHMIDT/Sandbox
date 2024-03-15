#include "CompuFluidDyna.hpp"


// Standard lib
#include <numbers>
#include <vector>

// Algo headers
#include "Math/Field.hpp"
#include "Math/Functions.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Incompressible Navier Stokes
// ∂vel/∂t = - (vel · ∇) vel + visco ∇²vel − 1/ρ ∇press + f
// ∇ · vel = 0
void CompuFluidDyna::RunSimulationStep() {
  // Get simulation parameters
  const int maxIter= std::max(D.UI[SolvMaxIter_____].I(), 0);
  const float timestep= std::max(D.UI[TimeStep________].F(), 0.0f);
  const float coeffDiffu= std::max(D.UI[CoeffDiffuS_____].F(), 0.0f);
  const float coeffVisco= std::max(D.UI[CoeffDiffuV_____].F(), 0.0f);
  const float coeffVorti= D.UI[CoeffVorti______].F();
  simTime+= D.UI[TimeStep________].F();

  // Update Solid boolean field based on continuous porosity value
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  UpdateSolidFieldFromPoros();
  if (D.UI[VerboseLevel____].I() >= 1) printf("UpdateSolid %f ", Timer::PopTimer());

  // Update periodic smoke in inlet
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ApplyBC(FieldID::IDSmok, Smok);
  if (D.UI[VerboseLevel____].I() >= 1) printf("ApplyBC %f ", Timer::PopTimer());

  // External forces
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ExternalGravityForce();
  ExternalDarcyForce();
  if (D.UI[VerboseLevel____].I() >= 1) printf("ExternF %f ", Timer::PopTimer());

  // Advection step with MacCormack backtracking
  if (D.UI[CoeffAdvec______].B()) {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    // Smoke advection
    // smo ⇐ smo - Δt (vel · ∇) smo
    AdvectField(FieldID::IDSmok, timestep, VelX, VelY, VelZ, Smok);
    // Non linear velocity self-advection
    // vel ⇐ vel - Δt (vel · ∇) vel
    std::vector<std::vector<std::vector<float>>> oldVelX= VelX;
    std::vector<std::vector<std::vector<float>>> oldVelY= VelY;
    std::vector<std::vector<std::vector<float>>> oldVelZ= VelZ;
    if (nX > 1) AdvectField(FieldID::IDVelX, timestep, oldVelX, oldVelY, oldVelZ, VelX);
    if (nY > 1) AdvectField(FieldID::IDVelY, timestep, oldVelX, oldVelY, oldVelZ, VelY);
    if (nZ > 1) AdvectField(FieldID::IDVelZ, timestep, oldVelX, oldVelY, oldVelZ, VelZ);
    if (D.UI[VerboseLevel____].I() >= 1) printf("Advect %f ", Timer::PopTimer());
  }

  // Smoke diffusion step
  // Implicit solve    (Id - diffu Δt ∇²) smo = smo
  if (D.UI[CoeffDiffuS_____].B()) {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    const std::vector<std::vector<std::vector<float>>> oldSmoke= Smok;
    LinearSolve(FieldID::IDSmok, maxIter, timestep, coeffDiffu, oldSmoke, Smok);
    if (D.UI[VerboseLevel____].I() >= 1) printf("DiffuS %f ", Timer::PopTimer());
  }

  // Velocity viscosity step
  // Implicit solve    (Id - visco Δt ∇²) vel = vel
  if (D.UI[CoeffDiffuV_____].B()) {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    const std::vector<std::vector<std::vector<float>>> oldVelX= VelX;
    const std::vector<std::vector<std::vector<float>>> oldVelY= VelY;
    const std::vector<std::vector<std::vector<float>>> oldVelZ= VelZ;
    if (nX > 1) LinearSolve(FieldID::IDVelX, maxIter, timestep, coeffVisco, oldVelX, VelX);
    if (nY > 1) LinearSolve(FieldID::IDVelY, maxIter, timestep, coeffVisco, oldVelY, VelY);
    if (nZ > 1) LinearSolve(FieldID::IDVelZ, maxIter, timestep, coeffVisco, oldVelZ, VelZ);
    if (D.UI[VerboseLevel____].I() >= 1) printf("DiffuV %f ", Timer::PopTimer());
  }

  // Vorticity step
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffVorti______].B()) VorticityConfinement(timestep, coeffVorti, VelX, VelY, VelZ);
  if (D.UI[VerboseLevel____].I() >= 1) printf("VortConfin %f ", Timer::PopTimer());

  // Projection step
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffProj_______].B()) ProjectField(maxIter, timestep, VelX, VelY, VelZ);
  if (D.UI[VerboseLevel____].I() >= 1) printf("Proj %f ", Timer::PopTimer());

  // Compute field data for display
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ComputeVectorFieldDivergence(VelX, VelY, VelZ, Dive);
  ComputeVectorFieldCurl(VelX, VelY, VelZ, CurX, CurY, CurZ);
  ComputeVectorFieldNorm(CurX, CurY, CurZ, Vort);
  if (D.UI[VerboseLevel____].I() >= 1) printf("DivCurlVort %f ", Timer::PopTimer());

  // Move particles for display
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ComputeParticlesMovement();
  if (D.UI[VerboseLevel____].I() >= 1) printf("MovPartic %f ", Timer::PopTimer());

  // Update the plots
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  UpdateUIData();
  if (D.UI[VerboseLevel____].I() >= 1) printf("UIData %f ", Timer::PopTimer());

  if (D.UI[VerboseLevel____].I() >= 1) printf("\n");
}


// Update Solid boolean field based on continuous porosity value
void CompuFluidDyna::UpdateSolidFieldFromPoros() {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        Solid[x][y][z]= (Poros[x][y][z] < D.UI[PorosMinThresh__].F());
}


// Apply boundary conditions enforcing fixed values to fields
void CompuFluidDyna::ApplyBC(const int iFieldID, std::vector<std::vector<std::vector<float>>>& ioField) {
  // Sweep through the field
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Set forced value
        if (Solid[x][y][z] && iFieldID == FieldID::IDSmok) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDVelX) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDVelY) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDVelZ) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDPres) ioField[x][y][z]= 0.0f;
        if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) ioField[x][y][z]= SmokForced[x][y][z] * std::cos(simTime * 2.0f * std::numbers::pi / D.UI[BCSmokTime______].F());
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelX) ioField[x][y][z]= VelXForced[x][y][z];
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelY) ioField[x][y][z]= VelYForced[x][y][z];
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelZ) ioField[x][y][z]= VelZForced[x][y][z];
        if (PreBC[x][y][z] && iFieldID == FieldID::IDPres) ioField[x][y][z]= PresForced[x][y][z];
      }
    }
  }
}


// Add external gravity/buoyancy forces to velocity field
// vel ⇐ vel + g * smo * Δt / ρ
void CompuFluidDyna::ExternalGravityForce() {
  if (D.UI[CoeffGravi______].F() != 0.0) {
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (Solid[x][y][z] || VelBC[x][y][z]) continue;
          VelZ[x][y][z]+= D.UI[CoeffGravi______].F() * Smok[x][y][z] * D.UI[TimeStep________].F() / fluidDensity;
        }
      }
    }
  }
}


// Add Darcy penalization forces from design shape to velocity field
// Implicit scheme: vel ⇐ vel / (1 + PenalInterpo(poros) * Δt / ρ)
void CompuFluidDyna::ExternalDarcyForce() {
  if (D.UI[DarcyMaxResist__].F() - D.UI[DarcyMinResist__].F() > 0.0f) {
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (Solid[x][y][z] || VelBC[x][y][z]) continue;
          const float penalVal= Functions::PenalInterpo(Poros[x][y][z], D.UI[DarcyPenal______].F(), false);
          const float resistVal= D.UI[DarcyMinResist__].F() + (D.UI[DarcyMaxResist__].F() - D.UI[DarcyMinResist__].F()) * (1.0f - penalVal);
          // Apply flow resistance force with implicit scheme
          // vel ⇐ vel / (1 + PenalInterpo(poros) * Δt / ρ)
          VelX[x][y][z]/= 1.0 + resistVal * D.UI[TimeStep________].F() / fluidDensity;
          VelY[x][y][z]/= 1.0 + resistVal * D.UI[TimeStep________].F() / fluidDensity;
          VelZ[x][y][z]/= 1.0 + resistVal * D.UI[TimeStep________].F() / fluidDensity;
          // // Apply flow resistance force with explicit scheme
          // // vel ⇐ vel - vel * PenalInterpo(poros) * Δt / ρ
          // VelX[x][y][z]-= VelX[x][y][z] * resistVal * D.UI[TimeStep________].F() / fluidDensity;
          // VelY[x][y][z]-= VelY[x][y][z] * resistVal * D.UI[TimeStep________].F() / fluidDensity;
          // VelZ[x][y][z]-= VelZ[x][y][z] * resistVal * D.UI[TimeStep________].F() / fluidDensity;
        }
      }
    }
  }
}


// Project velocity field into a solenoidal/divergence-free field
// Solve for pressure in pressure Poisson equation (negative Laplacian for positive diagonal in solve)
// (-∇²) press = -(ρ / Δt) ∇ · vel
// Update velocity field by subtracting gradient of pressure
// vel ⇐ vel - (Δt / ρ) × ∇ press
// References for pressure poisson equation and incompressiblity projection
// https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
// https://mycourses.aalto.fi/pluginfile.php/891524/mod_folder/content/0/Lecture03_Pressure.pdf
// https://barbagroup.github.io/essential_skills_RRC/numba/4/#application-pressure-poisson-equation
// http://www.thevisualroom.com/poisson_for_pressure.html
// https://github.com/barbagroup/CFDPython
void CompuFluidDyna::ProjectField(const int iMaxIter, const float iTimeStep,
                                  std::vector<std::vector<std::vector<float>>>& ioVelX,
                                  std::vector<std::vector<std::vector<float>>>& ioVelY,
                                  std::vector<std::vector<std::vector<float>>>& ioVelZ) {
  // Compute the RHS using the divergence of velocity
  std::vector<std::vector<std::vector<float>>> RHS= Field::AllocField3D(nX, nY, nZ, 0.0f);
  ComputeVectorFieldDivergence(VelX, VelY, VelZ, RHS);
  if (D.UI[TestPoroDiv_____].I() == 1) {
    for (int x= 0; x < nX; x++)
      for (int y= 0; y < nY; y++)
        for (int z= 0; z < nZ; z++)
          RHS[x][y][z]*= Poros[x][y][z];
  }
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        RHS[x][y][z]*= -(fluidDensity / D.UI[TimeStep________].F());
  ApplyBC(FieldID::IDPres, RHS);

  // Optionally reset pressure guess to test convergence
  if (D.UI[CoeffProj_______].I() == 2) {
    Pres= Field::AllocField3D(nX, nY, nZ, 0.0f);
    ApplyBC(FieldID::IDPres, Pres);
  }

  // Solve for pressure in the pressure Poisson equation
  LinearSolve(FieldID::IDPres, iMaxIter, iTimeStep, 0.0f, RHS, Pres);

  // Update velocities based on local pressure gradient
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z] || VelBC[x][y][z]) continue;
        // Subtract pressure gradient to remove divergence
        if (D.UI[TestPoroProj____].I() == 1) {
          if (x - 1 >= 0 && !Solid[x - 1][y][z]) ioVelX[x][y][z]-= Poros[x - 1][y][z] * iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x - 1][y][z]) / (2.0f * voxSize);
          if (y - 1 >= 0 && !Solid[x][y - 1][z]) ioVelY[x][y][z]-= Poros[x][y - 1][z] * iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x][y - 1][z]) / (2.0f * voxSize);
          if (z - 1 >= 0 && !Solid[x][y][z - 1]) ioVelZ[x][y][z]-= Poros[x][y][z - 1] * iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x][y][z - 1]) / (2.0f * voxSize);
          if (x + 1 < nX && !Solid[x + 1][y][z]) ioVelX[x][y][z]-= Poros[x + 1][y][z] * iTimeStep / fluidDensity * (Pres[x + 1][y][z] - Pres[x][y][z]) / (2.0f * voxSize);
          if (y + 1 < nY && !Solid[x][y + 1][z]) ioVelY[x][y][z]-= Poros[x][y + 1][z] * iTimeStep / fluidDensity * (Pres[x][y + 1][z] - Pres[x][y][z]) / (2.0f * voxSize);
          if (z + 1 < nZ && !Solid[x][y][z + 1]) ioVelZ[x][y][z]-= Poros[x][y][z + 1] * iTimeStep / fluidDensity * (Pres[x][y][z + 1] - Pres[x][y][z]) / (2.0f * voxSize);
        }
        else {
          if (x - 1 >= 0 && !Solid[x - 1][y][z]) ioVelX[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x - 1][y][z]) / (2.0f * voxSize);
          if (y - 1 >= 0 && !Solid[x][y - 1][z]) ioVelY[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x][y - 1][z]) / (2.0f * voxSize);
          if (z - 1 >= 0 && !Solid[x][y][z - 1]) ioVelZ[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x][y][z - 1]) / (2.0f * voxSize);
          if (x + 1 < nX && !Solid[x + 1][y][z]) ioVelX[x][y][z]-= iTimeStep / fluidDensity * (Pres[x + 1][y][z] - Pres[x][y][z]) / (2.0f * voxSize);
          if (y + 1 < nY && !Solid[x][y + 1][z]) ioVelY[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y + 1][z] - Pres[x][y][z]) / (2.0f * voxSize);
          if (z + 1 < nZ && !Solid[x][y][z + 1]) ioVelZ[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z + 1] - Pres[x][y][z]) / (2.0f * voxSize);
        }
      }
    }
  }
}


// Trilinearly interpolate the field value at the given position
float CompuFluidDyna::TrilinearInterpolation(const float iPosX, const float iPosY, const float iPosZ,
                                             const std::vector<std::vector<std::vector<float>>>& iFieldRef) {
  // Get floor and ceil voxel indices
  const int x0= std::min(std::max((int)std::floor(iPosX), 0), nX - 1);
  const int y0= std::min(std::max((int)std::floor(iPosY), 0), nY - 1);
  const int z0= std::min(std::max((int)std::floor(iPosZ), 0), nZ - 1);
  const int x1= std::min(std::max((int)std::ceil(iPosX), 0), nX - 1);
  const int y1= std::min(std::max((int)std::ceil(iPosY), 0), nY - 1);
  const int z1= std::min(std::max((int)std::ceil(iPosZ), 0), nZ - 1);
  // Get floor and ceil voxel weights
  const float xWeight1= iPosX - (float)x0;
  const float yWeight1= iPosY - (float)y0;
  const float zWeight1= iPosZ - (float)z0;
  const float xWeight0= 1.0f - xWeight1;
  const float yWeight0= 1.0f - yWeight1;
  const float zWeight0= 1.0f - zWeight1;
  // Compute the weighted sum
  return iFieldRef[x0][y0][z0] * (xWeight0 * yWeight0 * zWeight0) +
         iFieldRef[x0][y0][z1] * (xWeight0 * yWeight0 * zWeight1) +
         iFieldRef[x0][y1][z0] * (xWeight0 * yWeight1 * zWeight0) +
         iFieldRef[x0][y1][z1] * (xWeight0 * yWeight1 * zWeight1) +
         iFieldRef[x1][y0][z0] * (xWeight1 * yWeight0 * zWeight0) +
         iFieldRef[x1][y0][z1] * (xWeight1 * yWeight0 * zWeight1) +
         iFieldRef[x1][y1][z0] * (xWeight1 * yWeight1 * zWeight0) +
         iFieldRef[x1][y1][z1] * (xWeight1 * yWeight1 * zWeight1);
}


// Apply semi-Lagrangian advection along the velocity field
// vel ⇐ vel - Δt (vel · ∇) vel
// smo ⇐ smo - Δt (vel · ∇) smo
// References for MacCormack backtracking scheme
// https://commons.wikimedia.org/wiki/File:Backtracking_maccormack.png
// https://physbam.stanford.edu/~fedkiw/papers/stanford2006-09.pdf
// https://github.com/NiallHornFX/StableFluids3D-GL/blob/master/src/fluidsolver3d.cpp
void CompuFluidDyna::AdvectField(const int iFieldID, const float iTimeStep,
                                 const std::vector<std::vector<std::vector<float>>>& iVelX,
                                 const std::vector<std::vector<std::vector<float>>>& iVelY,
                                 const std::vector<std::vector<std::vector<float>>>& iVelZ,
                                 std::vector<std::vector<std::vector<float>>>& ioField) {
  // Adjust the source field to make solid voxels have a value dependant on their non-solid neighbors
  std::vector<std::vector<std::vector<float>>> sourceField= ioField;
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (!Solid[x][y][z]) continue;
        int count= 0;
        float sum= 0.0f;
        if (x - 1 >= 0 && !Solid[x - 1][y][z] && ++count) sum+= ioField[x - 1][y][z];
        if (y - 1 >= 0 && !Solid[x][y - 1][z] && ++count) sum+= ioField[x][y - 1][z];
        if (z - 1 >= 0 && !Solid[x][y][z - 1] && ++count) sum+= ioField[x][y][z - 1];
        if (x + 1 < nX && !Solid[x + 1][y][z] && ++count) sum+= ioField[x + 1][y][z];
        if (y + 1 < nY && !Solid[x][y + 1][z] && ++count) sum+= ioField[x][y + 1][z];
        if (z + 1 < nZ && !Solid[x][y][z + 1] && ++count) sum+= ioField[x][y][z + 1];
        if (iFieldID == FieldID::IDSmok) sourceField[x][y][z]= (count > 0) ? sum / (float)count : 0.0f;
        if (iFieldID == FieldID::IDVelX) sourceField[x][y][z]= (count > 0) ? -sum / (float)count : 0.0f;
        if (iFieldID == FieldID::IDVelY) sourceField[x][y][z]= (count > 0) ? -sum / (float)count : 0.0f;
        if (iFieldID == FieldID::IDVelZ) sourceField[x][y][z]= (count > 0) ? -sum / (float)count : 0.0f;
      }
    }
  }
  // Sweep through the field
#pragma omp parallel for collapse(3) if (D.UI[Multithread_____].B())
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        AdvX[x][y][z]= AdvY[x][y][z]= AdvZ[x][y][z]= 0.0f;
        // Skip solid or fixed values
        if (Solid[x][y][z]) continue;
        if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) continue;
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelX) continue;
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelY) continue;
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelZ) continue;
        // Find source position for active voxel using naive linear backtracking scheme
        const Vec::Vec3<float> posEnd((float)x, (float)y, (float)z);
        const Vec::Vec3<float> velEnd(iVelX[x][y][z], iVelY[x][y][z], iVelZ[x][y][z]);
        Vec::Vec3<float> posBeg= posEnd - iTimeStep * velEnd / voxSize;
        // Iterative source position correction with 2nd order MacCormack scheme
        for (int iter= 0; iter < D.UI[CoeffAdvec______].I() - 1; iter++) {
          const float velBegX= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelX);
          const float velBegY= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelY);
          const float velBegZ= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelZ);
          const Vec::Vec3<float> velBeg(velBegX, velBegY, velBegZ);
          const Vec::Vec3<float> vecErr= posEnd - (posBeg + iTimeStep * velBeg / voxSize);
          if (vecErr.normSquared() <= D.UI[CoeffAdvecTol___].F() * D.UI[CoeffAdvecTol___].F()) break;
          posBeg= posBeg + vecErr / 2.0f;
        }
        // Save source vector for display
        AdvX[x][y][z]= posBeg[0] - posEnd[0];
        AdvY[x][y][z]= posBeg[1] - posEnd[1];
        AdvZ[x][y][z]= posBeg[2] - posEnd[2];
        // Trilinear interpolation at source position
        ioField[x][y][z]= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], sourceField);
      }
    }
  }
}


// Counteract energy dissipation and introduce turbulent-like behavior by amplifying vorticity on small scales
// https://github.com/awesson/stable-fluids/tree/master
// https://github.com/woeishi/StableFluids/blob/master/StableFluid3d.cpp
// vel ⇐ vel + Δt * TODO write formula
void CompuFluidDyna::VorticityConfinement(const float iTimeStep, const float iVortiCoeff,
                                          std::vector<std::vector<std::vector<float>>>& ioVelX,
                                          std::vector<std::vector<std::vector<float>>>& ioVelY,
                                          std::vector<std::vector<std::vector<float>>>& ioVelZ) {
  // Compute curl and vorticity from the velocity field
  ComputeVectorFieldCurl(ioVelX, ioVelY, ioVelZ, CurX, CurY, CurZ);
  ComputeVectorFieldNorm(CurX, CurY, CurZ, Vort);
  // Amplify non-zero vorticity
  if (iVortiCoeff > 0.0f) {
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (Solid[x][y][z] || VelBC[x][y][z]) continue;
          // Gradient of vorticity with zero derivative at solid interface or domain boundary
          Vec::Vec3<float> vortGrad(0.0f, 0.0f, 0.0f);
          if (x - 1 >= 0 && !Solid[x - 1][y][z]) vortGrad[0]+= (Vort[x][y][z] - Vort[x - 1][y][z]) / (2.0f * voxSize);
          if (y - 1 >= 0 && !Solid[x][y - 1][z]) vortGrad[1]+= (Vort[x][y][z] - Vort[x][y - 1][z]) / (2.0f * voxSize);
          if (z - 1 >= 0 && !Solid[x][y][z - 1]) vortGrad[2]+= (Vort[x][y][z] - Vort[x][y][z - 1]) / (2.0f * voxSize);
          if (x + 1 < nX && !Solid[x + 1][y][z]) vortGrad[0]+= (Vort[x + 1][y][z] - Vort[x][y][z]) / (2.0f * voxSize);
          if (y + 1 < nY && !Solid[x][y + 1][z]) vortGrad[1]+= (Vort[x][y + 1][z] - Vort[x][y][z]) / (2.0f * voxSize);
          if (z + 1 < nZ && !Solid[x][y][z + 1]) vortGrad[2]+= (Vort[x][y][z + 1] - Vort[x][y][z]) / (2.0f * voxSize);
          // Amplification of small scale vorticity by following current curl
          if (vortGrad.norm() > 0.0f) {
            const float dVort_dx_scaled= iVortiCoeff * vortGrad[0] / vortGrad.norm();
            const float dVort_dy_scaled= iVortiCoeff * vortGrad[1] / vortGrad.norm();
            const float dVort_dz_scaled= iVortiCoeff * vortGrad[2] / vortGrad.norm();
            ioVelX[x][y][z]+= iTimeStep * (dVort_dy_scaled * CurZ[x][y][z] - dVort_dz_scaled * CurY[x][y][z]);
            ioVelY[x][y][z]+= iTimeStep * (dVort_dz_scaled * CurX[x][y][z] - dVort_dx_scaled * CurZ[x][y][z]);
            ioVelZ[x][y][z]+= iTimeStep * (dVort_dx_scaled * CurY[x][y][z] - dVort_dy_scaled * CurX[x][y][z]);
          }
        }
      }
    }
  }
}


// Move particles along velocity field
void CompuFluidDyna::ComputeParticlesMovement() {
  // Allocate particles if needed
  ParticlesPos.resize(D.UI[ParticCount_____].I());
  ParticlesAge.resize(D.UI[ParticCount_____].I(), -1.0f);
  // Reinitialize expired particles
  for (int k= 0; k < (int)ParticlesPos.size(); k++) {
    if (ParticlesAge[k] >= D.UI[ParticDuration__].F() ||
        ParticlesPos[k][0] < D.boxMin[0] || ParticlesPos[k][0] > D.boxMax[0] ||
        ParticlesPos[k][1] < D.boxMin[1] || ParticlesPos[k][1] > D.boxMax[1] ||
        ParticlesPos[k][2] < D.boxMin[2] || ParticlesPos[k][2] > D.boxMax[2]) {
      const float newPosX= (nX > 1) ? (Random::Val(D.boxMin[0], D.boxMax[0])) : (D.boxMax[0]);
      const float newPosY= (nY > 1) ? (Random::Val(D.boxMin[1], D.boxMax[1])) : (D.boxMax[1]);
      const float newPosZ= (nZ > 1) ? (Random::Val(D.boxMin[2], D.boxMax[2])) : (D.boxMax[2]);
      ParticlesPos[k].set(newPosX, newPosY, newPosZ);
      ParticlesAge[k]= Random::Val(0.0f, D.UI[ParticDuration__].F());
    }
  }
  // Update particles positions
  for (int k= 0; k < (int)ParticlesPos.size(); k++) {
    const float particleRelPosX= (float)nX * (ParticlesPos[k][0] - D.boxMin[0]) / (D.boxMax[0] - D.boxMin[0]);
    const float particleRelPosY= (float)nY * (ParticlesPos[k][1] - D.boxMin[1]) / (D.boxMax[1] - D.boxMin[1]);
    const float particleRelPosZ= (float)nZ * (ParticlesPos[k][2] - D.boxMin[2]) / (D.boxMax[2] - D.boxMin[2]);
    const int x= std::min(std::max(int(std::floor(particleRelPosX)), 0), nX - 1);
    const int y= std::min(std::max(int(std::floor(particleRelPosY)), 0), nY - 1);
    const int z= std::min(std::max(int(std::floor(particleRelPosZ)), 0), nZ - 1);
    ParticlesPos[k]+= 1.0f / 60.0f * Vec::Vec3<float>(VelX[x][y][z], VelY[x][y][z], VelZ[x][y][z]);
    ParticlesAge[k]+= 1.0f / 60.0f;
  }
}


// Compute the divergence of a vector field
// div= ∇ · vec
// Values in the solid are zero
// Values for the faces at the domain boundary are equal to the border cell (continuity)
// Values for the faces at the solid interface are equal to zero (mirror)
void CompuFluidDyna::ComputeVectorFieldDivergence(const std::vector<std::vector<std::vector<float>>>& iVecX,
                                                  const std::vector<std::vector<std::vector<float>>>& iVecY,
                                                  const std::vector<std::vector<std::vector<float>>>& iVecZ,
                                                  std::vector<std::vector<std::vector<float>>>& oDiv) {
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        oDiv[x][y][z]= 0.0f;
        if (Solid[x][y][z]) continue;
        // Classical linear interpolation for face values with same value at domain boundary and zero value at solid interface
        const float vecXN= (x - 1 >= 0) ? ((Solid[x - 1][y][z]) ? (0.0f) : ((iVecX[x][y][z] + iVecX[x - 1][y][z]) / 2.0f)) : (iVecX[x][y][z]);
        const float vecYN= (y - 1 >= 0) ? ((Solid[x][y - 1][z]) ? (0.0f) : ((iVecY[x][y][z] + iVecY[x][y - 1][z]) / 2.0f)) : (iVecY[x][y][z]);
        const float vecZN= (z - 1 >= 0) ? ((Solid[x][y][z - 1]) ? (0.0f) : ((iVecZ[x][y][z] + iVecZ[x][y][z - 1]) / 2.0f)) : (iVecZ[x][y][z]);
        const float vecXP= (x + 1 < nX) ? ((Solid[x + 1][y][z]) ? (0.0f) : ((iVecX[x + 1][y][z] + iVecX[x][y][z]) / 2.0f)) : (iVecX[x][y][z]);
        const float vecYP= (y + 1 < nY) ? ((Solid[x][y + 1][z]) ? (0.0f) : ((iVecY[x][y + 1][z] + iVecY[x][y][z]) / 2.0f)) : (iVecY[x][y][z]);
        const float vecZP= (z + 1 < nZ) ? ((Solid[x][y][z + 1]) ? (0.0f) : ((iVecZ[x][y][z + 1] + iVecZ[x][y][z]) / 2.0f)) : (iVecZ[x][y][z]);
        // Divergence based on face values
        oDiv[x][y][z]= ((vecXP - vecXN) + (vecYP - vecYN) + (vecZP - vecZN)) / voxSize;
      }
    }
  }
}


// Compute the curl of a vector field
// curl= ∇ ⨯ vec
void CompuFluidDyna::ComputeVectorFieldCurl(const std::vector<std::vector<std::vector<float>>>& iVecX,
                                            const std::vector<std::vector<std::vector<float>>>& iVecY,
                                            const std::vector<std::vector<std::vector<float>>>& iVecZ,
                                            std::vector<std::vector<std::vector<float>>>& oCurlX,
                                            std::vector<std::vector<std::vector<float>>>& oCurlY,
                                            std::vector<std::vector<std::vector<float>>>& oCurlZ) {
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        oCurlX[x][y][z]= oCurlY[x][y][z]= oCurlZ[x][y][z]= 0.0f;
        if (Solid[x][y][z]) continue;
        // Compute vector cross derivatives considering BC at interface with solid
        float dVecy_dx= 0.0f, dVecz_dx= 0.0f, dVecx_dy= 0.0f, dVecz_dy= 0.0f, dVecx_dz= 0.0f, dVecy_dz= 0.0f;
        if (x - 1 >= 0 && x + 1 < nX) dVecy_dx= ((Solid[x + 1][y][z] ? iVecY[x][y][z] : iVecY[x + 1][y][z]) - (Solid[x - 1][y][z] ? iVecY[x][y][z] : iVecY[x - 1][y][z])) / 2.0f;
        if (x - 1 >= 0 && x + 1 < nX) dVecz_dx= ((Solid[x + 1][y][z] ? iVecZ[x][y][z] : iVecZ[x + 1][y][z]) - (Solid[x - 1][y][z] ? iVecZ[x][y][z] : iVecZ[x - 1][y][z])) / 2.0f;
        if (y - 1 >= 0 && y + 1 < nY) dVecx_dy= ((Solid[x][y + 1][z] ? iVecX[x][y][z] : iVecX[x][y + 1][z]) - (Solid[x][y - 1][z] ? iVecX[x][y][z] : iVecX[x][y - 1][z])) / 2.0f;
        if (y - 1 >= 0 && y + 1 < nY) dVecz_dy= ((Solid[x][y + 1][z] ? iVecZ[x][y][z] : iVecZ[x][y + 1][z]) - (Solid[x][y - 1][z] ? iVecZ[x][y][z] : iVecZ[x][y - 1][z])) / 2.0f;
        if (z - 1 >= 0 && z + 1 < nZ) dVecx_dz= ((Solid[x][y][z + 1] ? iVecX[x][y][z] : iVecX[x][y][z + 1]) - (Solid[x][y][z - 1] ? iVecX[x][y][z] : iVecX[x][y][z - 1])) / 2.0f;
        if (z - 1 >= 0 && z + 1 < nZ) dVecy_dz= ((Solid[x][y][z + 1] ? iVecY[x][y][z] : iVecY[x][y][z + 1]) - (Solid[x][y][z - 1] ? iVecY[x][y][z] : iVecY[x][y][z - 1])) / 2.0f;
        // Deduce curl
        oCurlX[x][y][z]= dVecz_dy - dVecy_dz;
        oCurlY[x][y][z]= dVecx_dz - dVecz_dx;
        oCurlZ[x][y][z]= dVecy_dx - dVecx_dy;
      }
    }
  }
}


// Compute the norm of a vector field
// norm= ‖vec‖₂
void CompuFluidDyna::ComputeVectorFieldNorm(const std::vector<std::vector<std::vector<float>>>& iVecX,
                                            const std::vector<std::vector<std::vector<float>>>& iVecY,
                                            const std::vector<std::vector<std::vector<float>>>& iVecZ,
                                            std::vector<std::vector<std::vector<float>>>& oNorm) {
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        oNorm[x][y][z]= 0.0f;
        if (Solid[x][y][z]) continue;
        oNorm[x][y][z]= std::sqrt(iVecX[x][y][z] * iVecX[x][y][z] + iVecY[x][y][z] * iVecY[x][y][z] + iVecZ[x][y][z] * iVecZ[x][y][z]);
      }
    }
  }
}
