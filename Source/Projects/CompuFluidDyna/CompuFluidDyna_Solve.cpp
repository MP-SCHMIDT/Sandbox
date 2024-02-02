#include "CompuFluidDyna.hpp"


// Standard lib
#include <numbers>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void CompuFluidDyna::RunSimulationStep() {
  // Get simulation parameters
  const int maxIter= std::max(D.UI[SolvMaxIter_____].I(), 0);
  const float timestep= D.UI[TimeStep________].F();
  const float coeffDiffu= std::max(D.UI[CoeffDiffuS_____].F(), 0.0f);
  const float coeffVisco= std::max(D.UI[CoeffDiffuV_____].F(), 0.0f);
  const float coeffVorti= D.UI[CoeffVorti______].F();
  simTime+= D.UI[TimeStep________].F();

  // Precompute active voxel indices
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[PlotSolve_______].I() == 2 || D.UI[PlotSolve_______].I() == 3) {
    VoxIdx.clear();
    VoxIdx.resize(5);
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (Solid[x][y][z]) continue;
          if (!SmoBC[x][y][z]) VoxIdx[FieldID::IDSmok].push_back({x, y, z});
          if (!VelBC[x][y][z]) VoxIdx[FieldID::IDVelX].push_back({x, y, z});
          if (!VelBC[x][y][z]) VoxIdx[FieldID::IDVelY].push_back({x, y, z});
          if (!VelBC[x][y][z]) VoxIdx[FieldID::IDVelZ].push_back({x, y, z});
          if (!PreBC[x][y][z]) VoxIdx[FieldID::IDPres].push_back({x, y, z});
        }
      }
    }
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("ActiIdx %f ", Timer::PopTimer());

  // Update periodic smoke in inlet
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ApplyBC(FieldID::IDSmok, Smok);
  if (D.UI[VerboseLevel____].I() >= 1) printf("ApplyBC %f ", Timer::PopTimer());

  // Incompressible Navier Stokes
  // ∂vel/∂t = - (vel · ∇) vel + visco ∇²vel − 1/ρ ∇press + f
  // ∇ · vel = 0

  // Advection steps
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffAdvec______].B()) {
    AdvectField(FieldID::IDSmok, timestep, VelX, VelY, VelZ, Smok);
  }
  if (D.UI[CoeffAdvec______].B()) {
    std::vector<std::vector<std::vector<float>>> oldVelX= VelX;
    std::vector<std::vector<std::vector<float>>> oldVelY= VelY;
    std::vector<std::vector<std::vector<float>>> oldVelZ= VelZ;
    if (nX > 1) AdvectField(FieldID::IDVelX, timestep, oldVelX, oldVelY, oldVelZ, VelX);
    if (nY > 1) AdvectField(FieldID::IDVelY, timestep, oldVelX, oldVelY, oldVelZ, VelY);
    if (nZ > 1) AdvectField(FieldID::IDVelZ, timestep, oldVelX, oldVelY, oldVelZ, VelZ);
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("Advect %f ", Timer::PopTimer());

  // Diffusion steps
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffDiffuS_____].B()) {
    // (Id - diffu Δt ∇²) smo = smo
    std::vector<std::vector<std::vector<float>>> oldSmoke= Smok;
    if (D.UI[SolvType________].I() == 0) {
      GaussSeidelSolve(FieldID::IDSmok, maxIter, timestep, true, coeffDiffu, oldSmoke, Smok);
    }
    else if (D.UI[SolvType________].I() == 1) {
      GradientDescentSolve(FieldID::IDSmok, maxIter, timestep, true, coeffDiffu, oldSmoke, Smok);
    }
    else {
      ConjugateGradientSolve(FieldID::IDSmok, maxIter, timestep, true, coeffDiffu, oldSmoke, Smok);
    }
  }
  if (D.UI[CoeffDiffuV_____].B()) {
    // (Id - visco Δt ∇²) vel = vel
    std::vector<std::vector<std::vector<float>>> oldVelX= VelX;
    std::vector<std::vector<std::vector<float>>> oldVelY= VelY;
    std::vector<std::vector<std::vector<float>>> oldVelZ= VelZ;
    if (D.UI[SolvType________].I() == 0) {
      if (nX > 1) GaussSeidelSolve(FieldID::IDVelX, maxIter, timestep, true, coeffVisco, oldVelX, VelX);
      if (nY > 1) GaussSeidelSolve(FieldID::IDVelY, maxIter, timestep, true, coeffVisco, oldVelY, VelY);
      if (nZ > 1) GaussSeidelSolve(FieldID::IDVelZ, maxIter, timestep, true, coeffVisco, oldVelZ, VelZ);
    }
    else if (D.UI[SolvType________].I() == 1) {
      if (nX > 1) GradientDescentSolve(FieldID::IDVelX, maxIter, timestep, true, coeffVisco, oldVelX, VelX);
      if (nY > 1) GradientDescentSolve(FieldID::IDVelY, maxIter, timestep, true, coeffVisco, oldVelY, VelY);
      if (nZ > 1) GradientDescentSolve(FieldID::IDVelZ, maxIter, timestep, true, coeffVisco, oldVelZ, VelZ);
    }
    else {
      if (nX > 1) ConjugateGradientSolve(FieldID::IDVelX, maxIter, timestep, true, coeffVisco, oldVelX, VelX);
      if (nY > 1) ConjugateGradientSolve(FieldID::IDVelY, maxIter, timestep, true, coeffVisco, oldVelY, VelY);
      if (nZ > 1) ConjugateGradientSolve(FieldID::IDVelZ, maxIter, timestep, true, coeffVisco, oldVelZ, VelZ);
    }
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("Diffu %f ", Timer::PopTimer());

  // Vorticity step
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffVorti______].B()) {
    VorticityConfinement(timestep, coeffVorti, VelX, VelY, VelZ);
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("VortConfi %f ", Timer::PopTimer());

  // External forces
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffGravi______].B()) {
    ExternalForces();
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("ExternFor %f ", Timer::PopTimer());

  // Projection step
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  if (D.UI[CoeffProj_______].B()) {
    ProjectField(maxIter, timestep, VelX, VelY, VelZ);
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("Project %f ", Timer::PopTimer());

  // Compute field data for display
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ComputeVelocityDivergence();
  if (D.UI[VerboseLevel____].I() >= 1) printf("Div %f ", Timer::PopTimer());
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  ComputeVelocityCurlVorticity();
  if (D.UI[VerboseLevel____].I() >= 1) printf("CurlVorti %f ", Timer::PopTimer());

  // TODO Compute fluid density to check if constant as it should be in incompressible case

  // TODO Test heuristic optimization of solid regions
  // https://open-research-europe.ec.europa.eu/articles/3-156

  // TODO switch to continusous Solid field and implement Brinkman penalization ?
  // https://arxiv.org/pdf/2302.14156.pdf
  // https://link.springer.com/article/10.1007/s00158-023-03570-4
  // TODO Test with flow separation scenarios ?
  // TODO Test with diagonal vs axis aligned tube ?
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


// Addition of one field to an other
void CompuFluidDyna::ImplicitFieldAdd(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                                      const std::vector<std::vector<std::vector<float>>>& iFieldB,
                                      std::vector<std::vector<std::vector<float>>>& oField) {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        oField[x][y][z]= iFieldA[x][y][z] + iFieldB[x][y][z];
}


// Subtraction of one field to an other
void CompuFluidDyna::ImplicitFieldSub(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                                      const std::vector<std::vector<std::vector<float>>>& iFieldB,
                                      std::vector<std::vector<std::vector<float>>>& oField) {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        oField[x][y][z]= iFieldA[x][y][z] - iFieldB[x][y][z];
}


// Multiplication of field by scalar
void CompuFluidDyna::ImplicitFieldScale(const float iVal,
                                        const std::vector<std::vector<std::vector<float>>>& iField,
                                        std::vector<std::vector<std::vector<float>>>& oField) {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        oField[x][y][z]= iField[x][y][z] * iVal;
}


// Dot product between two fields
float CompuFluidDyna::ImplicitFieldDotProd(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                                           const std::vector<std::vector<std::vector<float>>>& iFieldB) {
  float val= 0.0f;
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        val+= iFieldA[x][y][z] * iFieldB[x][y][z];
  return val;
}


// Perform a matrix-vector multiplication without explicitly assembling the Laplacian matrix
void CompuFluidDyna::ImplicitFieldLaplacianMatMult(const int iFieldID, const float iTimeStep,
                                                   const bool iDiffuMode, const float iDiffuCoeff, const bool iPrecondMode,
                                                   const std::vector<std::vector<std::vector<float>>>& iField,
                                                   std::vector<std::vector<std::vector<float>>>& oField) {
  // Precompute value
  const float diffuVal= iDiffuCoeff * iTimeStep / (voxSize * voxSize);
  if (D.UI[PlotSolve_______].I() == 1) {
    // Sweep through the field
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          // Skip solid or fixed values
          if (Solid[x][y][z]) continue;
          if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) continue;
          if (VelBC[x][y][z] && (iFieldID == FieldID::IDVelX || iFieldID == FieldID::IDVelY || iFieldID == FieldID::IDVelZ)) continue;
          if (PreBC[x][y][z] && iFieldID == FieldID::IDPres) continue;
          // Get count and sum of valid neighbors
          const int count= (x > 0) + (y > 0) + (z > 0) + (x < nX - 1) + (y < nY - 1) + (z < nZ - 1);
          float sum= 0.0f;
          if (!iPrecondMode) {
            const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelX ? -iField[x][y][z] : 0.0f);
            const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelY ? -iField[x][y][z] : 0.0f);
            const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelZ ? -iField[x][y][z] : 0.0f);
            if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : iField[x - 1][y][z];
            if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : iField[x + 1][y][z];
            if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : iField[x][y - 1][z];
            if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : iField[x][y + 1][z];
            if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : iField[x][y][z - 1];
            if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : iField[x][y][z + 1];
          }
          // Apply linear expression
          if (iDiffuMode) {
            if (iPrecondMode)
              oField[x][y][z]= 1.0f / (1.0f + diffuVal * (float)count) * iField[x][y][z];            //               [   -D*dt/(h*h)]
            else                                                                                     // [-D*dt/(h*h)] [1+4*D*dt/(h*h)] [-D*dt/(h*h)]
              oField[x][y][z]= (1.0f + diffuVal * (float)count) * iField[x][y][z] - diffuVal * sum;  //               [   -D*dt/(h*h)]
          }
          else {
            if (iPrecondMode)
              oField[x][y][z]= ((voxSize * voxSize) / (float)count) * iField[x][y][z];        //            [-1/(h*h)]
            else                                                                              // [-1/(h*h)] [ 4/(h*h)] [-1/(h*h)]
              oField[x][y][z]= ((float)count * iField[x][y][z] - sum) / (voxSize * voxSize);  //            [-1/(h*h)]
          }
        }
      }
    }
  }
  if (D.UI[PlotSolve_______].I() == 2) {
    // Sweep through the field
    for (auto idx : VoxIdx[iFieldID]) {
      const int x= idx[0];
      const int y= idx[1];
      const int z= idx[2];
      // Get count and sum of valid neighbors
      const int count= (x > 0) + (y > 0) + (z > 0) + (x < nX - 1) + (y < nY - 1) + (z < nZ - 1);
      float sum= 0.0f;
      if (!iPrecondMode) {
        const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelX ? -iField[x][y][z] : 0.0f);
        const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelY ? -iField[x][y][z] : 0.0f);
        const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelZ ? -iField[x][y][z] : 0.0f);
        if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : iField[x - 1][y][z];
        if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : iField[x + 1][y][z];
        if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : iField[x][y - 1][z];
        if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : iField[x][y + 1][z];
        if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : iField[x][y][z - 1];
        if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : iField[x][y][z + 1];
      }
      // Apply linear expression
      if (iDiffuMode) {
        if (iPrecondMode)
          oField[x][y][z]= 1.0f / (1.0f + diffuVal * (float)count) * iField[x][y][z];            //               [   -D*dt/(h*h)]
        else                                                                                     // [-D*dt/(h*h)] [1+4*D*dt/(h*h)] [-D*dt/(h*h)]
          oField[x][y][z]= (1.0f + diffuVal * (float)count) * iField[x][y][z] - diffuVal * sum;  //               [   -D*dt/(h*h)]
      }
      else {
        if (iPrecondMode)
          oField[x][y][z]= ((voxSize * voxSize) / (float)count) * iField[x][y][z];        //            [-1/(h*h)]
        else                                                                              // [-1/(h*h)] [ 4/(h*h)] [-1/(h*h)]
          oField[x][y][z]= ((float)count * iField[x][y][z] - sum) / (voxSize * voxSize);  //            [-1/(h*h)]
      }
    }
  }
  if (D.UI[PlotSolve_______].I() == 3) {
// Sweep through the field
#pragma omp parallel for
    for (auto idx : VoxIdx[iFieldID]) {
      const int x= idx[0], y= idx[1], z= idx[2];
      // Get count and sum of valid neighbors
      const int count= (x > 0) + (y > 0) + (z > 0) + (x < nX - 1) + (y < nY - 1) + (z < nZ - 1);
      float sum= 0.0f;
      if (!iPrecondMode) {
        const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelX ? -iField[x][y][z] : 0.0f);
        const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelY ? -iField[x][y][z] : 0.0f);
        const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelZ ? -iField[x][y][z] : 0.0f);
        if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : iField[x - 1][y][z];
        if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : iField[x + 1][y][z];
        if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : iField[x][y - 1][z];
        if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : iField[x][y + 1][z];
        if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : iField[x][y][z - 1];
        if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : iField[x][y][z + 1];
      }
      // Apply linear expression
      if (iDiffuMode) {
        if (iPrecondMode)
          oField[x][y][z]= 1.0f / (1.0f + diffuVal * (float)count) * iField[x][y][z];            //               [   -D*dt/(h*h)]
        else                                                                                     // [-D*dt/(h*h)] [1+4*D*dt/(h*h)] [-D*dt/(h*h)]
          oField[x][y][z]= (1.0f + diffuVal * (float)count) * iField[x][y][z] - diffuVal * sum;  //               [   -D*dt/(h*h)]
      }
      else {
        if (iPrecondMode)
          oField[x][y][z]= ((voxSize * voxSize) / (float)count) * iField[x][y][z];        //            [-1/(h*h)]
        else                                                                              // [-1/(h*h)] [ 4/(h*h)] [-1/(h*h)]
          oField[x][y][z]= ((float)count * iField[x][y][z] - sum) / (voxSize * voxSize);  //            [-1/(h*h)]
      }
    }
  }
}


// Solve linear system with Conjugate Gradient approach
// References for linear solvers and particularily PCG
// https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// https://services.math.duke.edu/~holee/math361-2020/lectures/Conjugate_gradients.pdf
// https://www3.nd.edu/~zxu2/acms60212-40212-S12/final_project/Linear_solvers_GPU.pdf
// https://github.com/awesson/stable-fluids/tree/master
// https://en.wikipedia.org/wiki/Conjugate_gradient_method
void CompuFluidDyna::ConjugateGradientSolve(const int iFieldID, const int iMaxIter, const float iTimeStep,
                                            const bool iDiffuMode, const float iDiffuCoeff,
                                            const std::vector<std::vector<std::vector<float>>>& iField,
                                            std::vector<std::vector<std::vector<float>>>& ioField) {
  // Prepare convergence plot
  if (D.UI[PlotSolve_______].B()) {
    D.Plot.resize(5);
    D.Plot[FieldID::IDSmok].name= "Diffu S";
    D.Plot[FieldID::IDVelX].name= "Diffu VX";
    D.Plot[FieldID::IDVelY].name= "Diffu VY";
    D.Plot[FieldID::IDVelZ].name= "Diffu VZ";
    D.Plot[FieldID::IDPres].name= "Proj  P";
    D.Plot[iFieldID].val.clear();
    D.Plot[iFieldID].isLog= true;
  }
  // Allocate fields
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> dField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // Compute residual error magnitude    r = b - A x    errNew = r · r
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  const float normRHS= ImplicitFieldDotProd(iField, iField);
  float errNew= errBeg;
  dField= rField;
  // Error plot
  if (D.UI[VerboseLevel____].I() > 2) {
    if (iFieldID == FieldID::IDSmok) printf("CG Diffu S  [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelX) printf("CG Diffu VX [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelY) printf("CG Diffu VY [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelZ) printf("CG Diffu VZ [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDPres) printf("CG Proj  P  [%.2e] ", normRHS);
    printf("%.2e ", errNew);
  }
  if (D.UI[PlotSolve_______].B())
    D.Plot[iFieldID].val.push_back(errNew);
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= D.UI[SolvTolAbs______].F()) break;
    if (errNew / normRHS <= std::max(D.UI[SolvTolRhs______].F(), 0.0f)) break;
    if (errNew / errBeg <= std::max(D.UI[SolvTolRel______].F(), 0.0f)) break;
    // q = A d
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, dField, qField);
    // alpha = errNew / (d^T q)
    const float denom= ImplicitFieldDotProd(dField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // x = x + alpha d
    ImplicitFieldScale(alpha, dField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ioField= t1Field;
    ApplyBC(iFieldID, ioField);
    // r = r - alpha q
    ImplicitFieldScale(alpha, qField, t0Field);
    ImplicitFieldSub(rField, t0Field, t1Field);
    rField= t1Field;
    // errNew = r^T r
    const float errOld= errNew;
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[VerboseLevel____].I() > 2)
      printf("%.2e ", errNew);
    if (D.UI[PlotSolve_______].B())
      D.Plot[iFieldID].val.push_back(errNew);
    // d = r + (errNew / errOld) * d
    ImplicitFieldScale(errNew / errOld, dField, t0Field);
    ImplicitFieldAdd(rField, t0Field, dField);
  }
  if (D.UI[VerboseLevel____].I() > 2)
    printf("\n");
  // Error plot
  if (D.UI[PlotSolve_______].B()) {
    if (iFieldID == FieldID::IDSmok) Dum0= rField;
    if (iFieldID == FieldID::IDVelX) Dum1= rField;
    if (iFieldID == FieldID::IDVelY) Dum2= rField;
    if (iFieldID == FieldID::IDVelZ) Dum3= rField;
    if (iFieldID == FieldID::IDPres) Dum4= rField;
  }
}


// Solve linear system with Gradient descent approach
// Reference on page 55 of https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
void CompuFluidDyna::GradientDescentSolve(const int iFieldID, const int iMaxIter, const float iTimeStep,
                                          const bool iDiffuMode, const float iDiffuCoeff,
                                          const std::vector<std::vector<std::vector<float>>>& iField,
                                          std::vector<std::vector<std::vector<float>>>& ioField) {
  // Prepare convergence plot
  if (D.UI[PlotSolve_______].B()) {
    D.Plot.resize(5);
    D.Plot[FieldID::IDSmok].name= "Diffu S";
    D.Plot[FieldID::IDVelX].name= "Diffu VX";
    D.Plot[FieldID::IDVelY].name= "Diffu VY";
    D.Plot[FieldID::IDVelZ].name= "Diffu VZ";
    D.Plot[FieldID::IDPres].name= "Proj  P";
    D.Plot[iFieldID].val.clear();
    D.Plot[iFieldID].isLog= true;
  }
  // Allocate fields
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // Compute residual error magnitude    r = b - A x    errNew = r · r
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  const float normRHS= ImplicitFieldDotProd(iField, iField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[VerboseLevel____].I() > 2) {
    if (iFieldID == FieldID::IDSmok) printf("GD Diffu S  [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelX) printf("GD Diffu VX [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelY) printf("GD Diffu VY [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelZ) printf("GD Diffu VZ [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDPres) printf("GD Proj  P  [%.2e] ", normRHS);
    printf("%.2e ", errNew);
  }
  if (D.UI[PlotSolve_______].B())
    D.Plot[iFieldID].val.push_back(errNew);
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= 0.0f) break;
    if (errNew <= D.UI[SolvTolAbs______].F()) break;
    if (errNew / normRHS <= std::max(D.UI[SolvTolRhs______].F(), 0.0f)) break;
    if (errNew / errBeg <= std::max(D.UI[SolvTolRel______].F(), 0.0f)) break;
    // q = A r
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, rField, qField);
    // alpha = errNew / (r^T q)
    const float denom= ImplicitFieldDotProd(rField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // x = x + alpha r
    ImplicitFieldScale(alpha, rField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ioField= t1Field;
    ApplyBC(iFieldID, ioField);
    // r = r - alpha q
    ImplicitFieldScale(alpha, qField, t0Field);
    ImplicitFieldSub(rField, t0Field, t1Field);
    rField= t1Field;
    // errNew = r^T r
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[VerboseLevel____].I() > 2)
      printf("%.2e ", errNew);
    if (D.UI[PlotSolve_______].B())
      D.Plot[iFieldID].val.push_back(errNew);
  }
  if (D.UI[VerboseLevel____].I() > 2)
    printf("\n");
  // Error plot
  if (D.UI[PlotSolve_______].B()) {
    if (iFieldID == FieldID::IDSmok) Dum0= rField;
    if (iFieldID == FieldID::IDVelX) Dum1= rField;
    if (iFieldID == FieldID::IDVelY) Dum2= rField;
    if (iFieldID == FieldID::IDVelZ) Dum3= rField;
    if (iFieldID == FieldID::IDPres) Dum4= rField;
  }
}


// Solve linear system with an iterative Gauss Seidel scheme
// Solving each equation sequentially by cascading latest solution to next equation
// Run a forward and backward pass in parallel to avoid element ordering bias
// Apply successive overrelaxation coefficient to accelerate convergence
void CompuFluidDyna::GaussSeidelSolve(const int iFieldID, const int iMaxIter, const float iTimeStep,
                                      const bool iDiffuMode, const float iDiffuCoeff,
                                      const std::vector<std::vector<std::vector<float>>>& iField,
                                      std::vector<std::vector<std::vector<float>>>& ioField) {
  // Prepare convergence plot
  if (D.UI[PlotSolve_______].B()) {
    D.Plot.resize(5);
    D.Plot[FieldID::IDSmok].name= "Diffu S";
    D.Plot[FieldID::IDVelX].name= "Diffu VX";
    D.Plot[FieldID::IDVelY].name= "Diffu VY";
    D.Plot[FieldID::IDVelZ].name= "Diffu VZ";
    D.Plot[FieldID::IDPres].name= "Proj  P";
    D.Plot[iFieldID].val.clear();
    D.Plot[iFieldID].isLog= true;
  }
  // Allocate fields
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<std::vector<float>>>> FieldT(2);
  // Compute residual error magnitude    r = b - A x    errNew = r · r
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  const float normRHS= ImplicitFieldDotProd(iField, iField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[VerboseLevel____].I() > 2) {
    if (iFieldID == FieldID::IDSmok) printf("GS Diffu S  [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelX) printf("GS Diffu VX [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelY) printf("GS Diffu VY [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelZ) printf("GS Diffu VZ [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDPres) printf("GS Proj  P  [%.2e] ", normRHS);
    printf("%.2e ", errNew);
  }
  if (D.UI[PlotSolve_______].B())
    D.Plot[iFieldID].val.push_back(errNew);
  // Precompute values
  const float diffuVal= iDiffuCoeff * iTimeStep / (voxSize * voxSize);
  const float coeffOverrelax= std::max(D.UI[SolvSOR_________].F(), 0.0f);
  // Iterate to solve with Gauss-Seidel scheme
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= D.UI[SolvTolAbs______].F()) break;
    if (errNew / normRHS <= std::max(D.UI[SolvTolRhs______].F(), 0.0f)) break;
    if (errNew / errBeg <= std::max(D.UI[SolvTolRel______].F(), 0.0f)) break;
    // Initialize fields for forward and backward passes
    FieldT[0]= ioField;
    FieldT[1]= ioField;
    // Execute the two passes in parallel
#pragma omp parallel for
    for (int k= 0; k < 2; k++) {
      // Set the loop settings for the current pass
      const int xBeg= (k == 0) ? 0 : nX - 1;
      const int yBeg= (k == 0) ? 0 : nY - 1;
      const int zBeg= (k == 0) ? 0 : nZ - 1;
      const int xEnd= (k == 0) ? nX : -1;
      const int yEnd= (k == 0) ? nY : -1;
      const int zEnd= (k == 0) ? nZ : -1;
      const int xInc= (k == 0) ? 1 : -1;
      const int yInc= (k == 0) ? 1 : -1;
      const int zInc= (k == 0) ? 1 : -1;
      // Sweep through the field
      for (int x= xBeg; x != xEnd; x+= xInc) {
        for (int y= yBeg; y != yEnd; y+= yInc) {
          for (int z= zBeg; z != zEnd; z+= zInc) {
            // Skip solid or fixed values
            if (Solid[x][y][z]) continue;
            if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) continue;
            if (VelBC[x][y][z] && (iFieldID == FieldID::IDVelX || iFieldID == FieldID::IDVelY || iFieldID == FieldID::IDVelZ)) continue;
            if (PreBC[x][y][z] && iFieldID == FieldID::IDPres) continue;
            // Get count and sum of valid neighbors
            const int count= (x > 0) + (y > 0) + (z > 0) + (x < nX - 1) + (y < nY - 1) + (z < nZ - 1);
            float sum= 0.0f;
            const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (FieldT[k][x][y][z]) : (iFieldID == FieldID::IDVelX ? -FieldT[k][x][y][z] : 0.0f);
            const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (FieldT[k][x][y][z]) : (iFieldID == FieldID::IDVelY ? -FieldT[k][x][y][z] : 0.0f);
            const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (FieldT[k][x][y][z]) : (iFieldID == FieldID::IDVelZ ? -FieldT[k][x][y][z] : 0.0f);
            if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : FieldT[k][x - 1][y][z];
            if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : FieldT[k][x + 1][y][z];
            if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : FieldT[k][x][y - 1][z];
            if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : FieldT[k][x][y + 1][z];
            if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : FieldT[k][x][y][z - 1];
            if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : FieldT[k][x][y][z + 1];
            // Set new value according to coefficients and flags
            if (count > 0) {
              const float prevVal= FieldT[k][x][y][z];
              if (iDiffuMode) FieldT[k][x][y][z]= (iField[x][y][z] + diffuVal * sum) / (1.0f + diffuVal * (float)count);
              else FieldT[k][x][y][z]= ((voxSize * voxSize) * iField[x][y][z] + sum) / (float)count;
              FieldT[k][x][y][z]= prevVal + coeffOverrelax * (FieldT[k][x][y][z] - prevVal);
            }
          }
        }
      }
    }
    // Recombine forward and backward passes
    for (int x= 0; x < nX; x++)
      for (int y= 0; y < nY; y++)
        for (int z= 0; z < nZ; z++)
          ioField[x][y][z]= (FieldT[0][x][y][z] + FieldT[1][x][y][z]) / 2.0f;
    // Compute residual error magnitude    r = b - A x    errNew = r · r
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
    ApplyBC(iFieldID, t0Field);
    ImplicitFieldSub(iField, t0Field, rField);
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[VerboseLevel____].I() > 2)
      printf("%.2e ", errNew);
    if (D.UI[PlotSolve_______].B())
      D.Plot[iFieldID].val.push_back(errNew);
  }
  if (D.UI[VerboseLevel____].I() > 2)
    printf("\n");
  // Error plot
  if (D.UI[PlotSolve_______].B()) {
    if (iFieldID == FieldID::IDSmok) Dum0= rField;
    if (iFieldID == FieldID::IDVelX) Dum1= rField;
    if (iFieldID == FieldID::IDVelY) Dum2= rField;
    if (iFieldID == FieldID::IDVelZ) Dum3= rField;
    if (iFieldID == FieldID::IDPres) Dum4= rField;
  }
}


// Add external forces to velocity field
// vel ⇐ vel + Δt * F / ρ
void CompuFluidDyna::ExternalForces() {
  // Update velocities based on applied external forces
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z] || VelBC[x][y][z]) continue;
        VelZ[x][y][z]+= D.UI[TimeStep________].F() * D.UI[CoeffGravi______].F() * Smok[x][y][z] / fluidDensity;
      }
    }
  }
}


// Project velocity field into a solenoidal/divergence-free field
// 1. Compute RHS based on divergence
// RHS = -(ρ / Δt) × ∇ · vel
// 2. Solve for pressure in pressure Poisson equation
// (-∇²) press = RHS
// 3. Update velocity field by subtracting gradient of pressure
// vel ⇐ vel - (Δt / ρ) × ∇ press
// References for pressure poisson equation and incompressiblity projection
// https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
// https://mycourses.aalto.fi/pluginfile.php/891524/mod_folder/content/0/Lecture03_Pressure.pdf
// https://barbagroup.github.io/essential_skills_RRC/numba/4/#application-pressure-poisson-equation
// http://www.thevisualroom.com/poisson_for_pressure.html
// https://github.com/barbagroup/CFDPython
void CompuFluidDyna::ProjectField(const int iIter, const float iTimeStep,
                                  std::vector<std::vector<std::vector<float>>>& ioVelX,
                                  std::vector<std::vector<std::vector<float>>>& ioVelY,
                                  std::vector<std::vector<std::vector<float>>>& ioVelZ) {
  // Compute divergence for RHS
  ComputeVelocityDivergence();
  // Reset pressure guess to test convergence
  if (D.UI[CoeffProj_______].I() == 2) {
    Pres= Field::AllocField3D(nX, nY, nZ, 0.0f);
    ApplyBC(FieldID::IDPres, Pres);
  }
  // Solve for pressure in the pressure Poisson equation
  if (D.UI[SolvType________].I() == 0) {
    GaussSeidelSolve(FieldID::IDPres, iIter, iTimeStep, false, 0.0f, Dive, Pres);
  }
  else if (D.UI[SolvType________].I() == 1) {
    GradientDescentSolve(FieldID::IDPres, iIter, iTimeStep, false, 0.0f, Dive, Pres);
  }
  else {
    ConjugateGradientSolve(FieldID::IDPres, iIter, iTimeStep, false, 0.0f, Dive, Pres);
  }

  // Update velocities based on local pressure gradient
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z] || VelBC[x][y][z]) continue;
        // Subtract pressure gradient to remove divergence
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
  for (int x= 0; x < nX; x++) {
#pragma omp parallel for
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
        int correcMaxIter= std::max(D.UI[CoeffAdvec______].I() - 1, 0);
        for (int iter= 0; iter < correcMaxIter; iter++) {
          const float velBegX= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelX);
          const float velBegY= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelY);
          const float velBegZ= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelZ);
          const Vec::Vec3<float> velBeg(velBegX, velBegY, velBegZ);
          const Vec::Vec3<float> vecErr= posEnd - (posBeg + iTimeStep * velBeg / voxSize);
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
  ComputeVelocityCurlVorticity();
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


// Compute RHS of pressure poisson equation as negative divergence scaled by density and timestep
// https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
// RHS = -(ρ / Δt) × ∇ · vel
// TODO implement correction to avoid checkerboard due to odd-even decoupling in pure pressure driven flows
// References for Rhie Chow correction
// https://youtu.be/yqZ59Xn_aF8 Checkerboard oscillations
// https://youtu.be/PmEUiUB8ETk Deriving the correction
// https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2007/rhiechow.pdf OpenFOAM variant
// https://mustafabhotvawala.com/wp-content/uploads/2020/11/MB_rhieChow-1.pdf
void CompuFluidDyna::ComputeVelocityDivergence() {
  // // Precompute pressure gradient for Rhie and Chow correction
  // std::vector<std::vector<std::vector<float>>> PresGradX;
  // std::vector<std::vector<std::vector<float>>> PresGradY;
  // std::vector<std::vector<std::vector<float>>> PresGradZ;
  // if (iUseRhieChow) {
  //   PresGradX= Field::AllocField3D(nX, nY, nZ, 0.0f);
  //   PresGradY= Field::AllocField3D(nX, nY, nZ, 0.0f);
  //   PresGradZ= Field::AllocField3D(nX, nY, nZ, 0.0f);
  //   for (int x= 0; x < nX; x++) {
  //     for (int y= 0; y < nY; y++) {
  //       for (int z= 0; z < nZ; z++) {
  //         if (Solid[x][y][z]) continue;
  //         // Pressure gradient with zero derivative at solid interface or domain boundary
  //         if (x - 1 >= 0 && !Solid[x - 1][y][z]) PresGradX[x][y][z]+= (Pres[x][y][z] - Pres[x - 1][y][z]) / (2.0f * voxSize);
  //         if (y - 1 >= 0 && !Solid[x][y - 1][z]) PresGradY[x][y][z]+= (Pres[x][y][z] - Pres[x][y - 1][z]) / (2.0f * voxSize);
  //         if (z - 1 >= 0 && !Solid[x][y][z - 1]) PresGradZ[x][y][z]+= (Pres[x][y][z] - Pres[x][y][z - 1]) / (2.0f * voxSize);
  //         if (x + 1 < nX && !Solid[x + 1][y][z]) PresGradX[x][y][z]+= (Pres[x + 1][y][z] - Pres[x][y][z]) / (2.0f * voxSize);
  //         if (y + 1 < nY && !Solid[x][y + 1][z]) PresGradY[x][y][z]+= (Pres[x][y + 1][z] - Pres[x][y][z]) / (2.0f * voxSize);
  //         if (z + 1 < nZ && !Solid[x][y][z + 1]) PresGradZ[x][y][z]+= (Pres[x][y][z + 1] - Pres[x][y][z]) / (2.0f * voxSize);
  //       }
  //     }
  //   }
  // }
  // Compute divergence of velocity field
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z]) Dive[x][y][z]= 0.0f;
        if (PreBC[x][y][z]) Dive[x][y][z]= PresForced[x][y][z];
        if (Solid[x][y][z] || PreBC[x][y][z]) continue;
        // Classical linear interpolation for face velocities with same velocity at domain boundary and zero velocity at solid interface
        float velXN= (x - 1 >= 0) ? ((Solid[x - 1][y][z]) ? (0.0f) : ((VelX[x][y][z] + VelX[x - 1][y][z]) / 2.0f)) : (VelX[x][y][z]);
        float velYN= (y - 1 >= 0) ? ((Solid[x][y - 1][z]) ? (0.0f) : ((VelY[x][y][z] + VelY[x][y - 1][z]) / 2.0f)) : (VelY[x][y][z]);
        float velZN= (z - 1 >= 0) ? ((Solid[x][y][z - 1]) ? (0.0f) : ((VelZ[x][y][z] + VelZ[x][y][z - 1]) / 2.0f)) : (VelZ[x][y][z]);
        float velXP= (x + 1 < nX) ? ((Solid[x + 1][y][z]) ? (0.0f) : ((VelX[x + 1][y][z] + VelX[x][y][z]) / 2.0f)) : (VelX[x][y][z]);
        float velYP= (y + 1 < nY) ? ((Solid[x][y + 1][z]) ? (0.0f) : ((VelY[x][y + 1][z] + VelY[x][y][z]) / 2.0f)) : (VelY[x][y][z]);
        float velZP= (z + 1 < nZ) ? ((Solid[x][y][z + 1]) ? (0.0f) : ((VelZ[x][y][z + 1] + VelZ[x][y][z]) / 2.0f)) : (VelZ[x][y][z]);
        // // Rhie and Chow correction terms
        // if (iUseRhieChow) {
        //   // Subtract pressure gradients with neighboring cells
        //   velXN-= D.UI[CoeffProj1______].F() * ((x - 1 >= 0 && !Solid[x - 1][y][z]) ? ((Pres[x][y][z] - Pres[x - 1][y][z]) / voxSize) : (0.0f));
        //   velYN-= D.UI[CoeffProj1______].F() * ((y - 1 >= 0 && !Solid[x][y - 1][z]) ? ((Pres[x][y][z] - Pres[x][y - 1][z]) / voxSize) : (0.0f));
        //   velZN-= D.UI[CoeffProj1______].F() * ((z - 1 >= 0 && !Solid[x][y][z - 1]) ? ((Pres[x][y][z] - Pres[x][y][z - 1]) / voxSize) : (0.0f));
        //   velXP-= D.UI[CoeffProj1______].F() * ((x + 1 < nX && !Solid[x + 1][y][z]) ? ((Pres[x + 1][y][z] - Pres[x][y][z]) / voxSize) : (0.0f));
        //   velYP-= D.UI[CoeffProj1______].F() * ((y + 1 < nY && !Solid[x][y + 1][z]) ? ((Pres[x][y + 1][z] - Pres[x][y][z]) / voxSize) : (0.0f));
        //   velZP-= D.UI[CoeffProj1______].F() * ((z + 1 < nZ && !Solid[x][y][z + 1]) ? ((Pres[x][y][z + 1] - Pres[x][y][z]) / voxSize) : (0.0f));
        //   // Add Linear interpolations of pressure gradients with neighboring cells
        //   velXN+= D.UI[CoeffProj2______].F() * ((x - 1 >= 0) ? ((PresGradX[x][y][z] + PresGradX[x - 1][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velYN+= D.UI[CoeffProj2______].F() * ((y - 1 >= 0) ? ((PresGradY[x][y][z] + PresGradY[x][y - 1][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velZN+= D.UI[CoeffProj2______].F() * ((z - 1 >= 0) ? ((PresGradZ[x][y][z] + PresGradZ[x][y][z - 1]) / 2.0f) : (PresGradX[x][y][z]));
        //   velXP+= D.UI[CoeffProj2______].F() * ((x + 1 < nX) ? ((PresGradX[x + 1][y][z] + PresGradX[x][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velYP+= D.UI[CoeffProj2______].F() * ((y + 1 < nY) ? ((PresGradY[x][y + 1][z] + PresGradY[x][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velZP+= D.UI[CoeffProj2______].F() * ((z + 1 < nZ) ? ((PresGradZ[x][y][z + 1] + PresGradZ[x][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        // }
        // Divergence based on face velocities scaled by density and timestep  (negated RHS and linear system to have positive diag coeffs)
        Dive[x][y][z]= -fluidDensity / D.UI[TimeStep________].F() * ((velXP - velXN) + (velYP - velYN) + (velZP - velZN)) / voxSize;
      }
    }
  }
}


// Compute curl and vorticity of current velocity field
// curl= ∇ ⨯ vel
// vort= ‖curl‖₂
void CompuFluidDyna::ComputeVelocityCurlVorticity() {
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        CurX[x][y][z]= CurY[x][y][z]= CurZ[x][y][z]= Vort[x][y][z]= 0.0f;
        if (Solid[x][y][z]) continue;
        // Compute velocity cross derivatives considering BC at interface with solid
        float dVely_dx= 0.0f, dVelz_dx= 0.0f, dVelx_dy= 0.0f, dVelz_dy= 0.0f, dVelx_dz= 0.0f, dVely_dz= 0.0f;
        if (x - 1 >= 0 && x + 1 < nX) dVely_dx= ((Solid[x + 1][y][z] ? VelY[x][y][z] : VelY[x + 1][y][z]) - (Solid[x - 1][y][z] ? VelY[x][y][z] : VelY[x - 1][y][z])) / 2.0f;
        if (x - 1 >= 0 && x + 1 < nX) dVelz_dx= ((Solid[x + 1][y][z] ? VelZ[x][y][z] : VelZ[x + 1][y][z]) - (Solid[x - 1][y][z] ? VelZ[x][y][z] : VelZ[x - 1][y][z])) / 2.0f;
        if (y - 1 >= 0 && y + 1 < nY) dVelx_dy= ((Solid[x][y + 1][z] ? VelX[x][y][z] : VelX[x][y + 1][z]) - (Solid[x][y - 1][z] ? VelX[x][y][z] : VelX[x][y - 1][z])) / 2.0f;
        if (y - 1 >= 0 && y + 1 < nY) dVelz_dy= ((Solid[x][y + 1][z] ? VelZ[x][y][z] : VelZ[x][y + 1][z]) - (Solid[x][y - 1][z] ? VelZ[x][y][z] : VelZ[x][y - 1][z])) / 2.0f;
        if (z - 1 >= 0 && z + 1 < nZ) dVelx_dz= ((Solid[x][y][z + 1] ? VelX[x][y][z] : VelX[x][y][z + 1]) - (Solid[x][y][z - 1] ? VelX[x][y][z] : VelX[x][y][z - 1])) / 2.0f;
        if (z - 1 >= 0 && z + 1 < nZ) dVely_dz= ((Solid[x][y][z + 1] ? VelY[x][y][z] : VelY[x][y][z + 1]) - (Solid[x][y][z - 1] ? VelY[x][y][z] : VelY[x][y][z - 1])) / 2.0f;
        // Deduce curl and vorticity
        CurX[x][y][z]= dVelz_dy - dVely_dz;
        CurY[x][y][z]= dVelx_dz - dVelz_dx;
        CurZ[x][y][z]= dVely_dx - dVelx_dy;
        Vort[x][y][z]= std::sqrt(CurX[x][y][z] * CurX[x][y][z] + CurY[x][y][z] * CurY[x][y][z] + CurZ[x][y][z] * CurZ[x][y][z]);
      }
    }
  }
}
