#include "CompuFluidDyna.hpp"


// Standard lib
#include <vector>

// Algo headers
#include "Math/Field.hpp"
#include "Math/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


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
void CompuFluidDyna::ImplicitFieldLaplacianMatMult(const int iFieldID, const float iTimeStep, const float iDiffuCoeff,
                                                   const std::vector<std::vector<std::vector<float>>>& iField,
                                                   std::vector<std::vector<std::vector<float>>>& oField) {
  // Precompute values
  const float voxSizeSqr= voxSize * voxSize;
  const float diffuVal= iDiffuCoeff * iTimeStep / voxSizeSqr;
  // Sweep through the field
#pragma omp parallel for collapse(3) if (D.UI[Multithread_____].B())
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
        const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelX ? -iField[x][y][z] : 0.0f);
        const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelY ? -iField[x][y][z] : 0.0f);
        const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelZ ? -iField[x][y][z] : 0.0f);
        if (D.UI[TestPoroLapla___].B()) {
          if (x - 1 >= 0) sum+= (1.0f - Poros[x - 1][y][z]) * xBCVal + Poros[x - 1][y][z] * iField[x - 1][y][z];
          if (x + 1 < nX) sum+= (1.0f - Poros[x + 1][y][z]) * xBCVal + Poros[x + 1][y][z] * iField[x + 1][y][z];
          if (y - 1 >= 0) sum+= (1.0f - Poros[x][y - 1][z]) * yBCVal + Poros[x][y - 1][z] * iField[x][y - 1][z];
          if (y + 1 < nY) sum+= (1.0f - Poros[x][y + 1][z]) * yBCVal + Poros[x][y + 1][z] * iField[x][y + 1][z];
          if (z - 1 >= 0) sum+= (1.0f - Poros[x][y][z - 1]) * zBCVal + Poros[x][y][z - 1] * iField[x][y][z - 1];
          if (z + 1 < nZ) sum+= (1.0f - Poros[x][y][z + 1]) * zBCVal + Poros[x][y][z + 1] * iField[x][y][z + 1];
        }
        else {
          if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : iField[x - 1][y][z];
          if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : iField[x + 1][y][z];
          if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : iField[x][y - 1][z];
          if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : iField[x][y + 1][z];
          if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : iField[x][y][z - 1];
          if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : iField[x][y][z + 1];
        }
        // Apply linear expression
        if (iFieldID == FieldID::IDPres)
          oField[x][y][z]= ((float)count * iField[x][y][z] - sum) / voxSizeSqr;  // [-1/(h*h)]  [ 2/(h*h)]  [-1/(h*h)]
        else
          oField[x][y][z]= (1.0f + (float)count * diffuVal) * iField[x][y][z] - diffuVal * sum;  // [-D*dt/(h*h)]  [1+2*D*dt/(h*h)]  [-D*dt/(h*h)]
      }
    }
  }
}


// Perform a matrix-vector multiplication without explicitly assembling the Laplacian matrix
void CompuFluidDyna::ImplicitFieldLaplacianDiagPrecond(const int iFieldID, const float iTimeStep, const float iDiffuCoeff,
                                                       const std::vector<std::vector<std::vector<float>>>& iField,
                                                       std::vector<std::vector<std::vector<float>>>& oField) {
  // Precompute values
  const float voxSizeSqr= voxSize * voxSize;
  const float diffuVal= iDiffuCoeff * iTimeStep / voxSizeSqr;
  // Sweep through the field
#pragma omp parallel for collapse(3) if (D.UI[Multithread_____].B())
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
        // Apply linear expression
        if (iFieldID == FieldID::IDPres)
          oField[x][y][z]= (voxSizeSqr / (float)count) * iField[x][y][z];
        else
          oField[x][y][z]= 1.0f / (1.0f + (float)count * diffuVal) * iField[x][y][z];
      }
    }
  }
}


// Wrapper function to call Gauss Seidel, Gradient Descent or Preconditioned Conjugate Gradient solver
void CompuFluidDyna::LinearSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                                 const std::vector<std::vector<std::vector<float>>>& iField,
                                 std::vector<std::vector<std::vector<float>>>& ioField) {
  // Initialize convergence plot
  if (D.UI[PlotSolve_______].B()) {
    D.Plot.push_back(PlotUI());
    if (iFieldID == FieldID::IDSmok) D.Plot[D.Plot.size() - 1].name= "Diffu S ";
    if (iFieldID == FieldID::IDVelX) D.Plot[D.Plot.size() - 1].name= "Diffu VX";
    if (iFieldID == FieldID::IDVelY) D.Plot[D.Plot.size() - 1].name= "Diffu VY";
    if (iFieldID == FieldID::IDVelZ) D.Plot[D.Plot.size() - 1].name= "Diffu VZ";
    if (iFieldID == FieldID::IDPres) D.Plot[D.Plot.size() - 1].name= "Proj  P ";
  }
  // Run the solver
  if (D.UI[SolvType________].I() == 0)
    GaussSeidelSolve(iFieldID, iMaxIter, iTimeStep, iDiffuCoeff, iField, ioField);
  else if (D.UI[SolvType________].I() == 1)
    GradientDescentSolve(iFieldID, iMaxIter, iTimeStep, iDiffuCoeff, iField, ioField);
  else if (D.UI[SolvType________].I() == 2)
    ConjugateGradientSolve(iFieldID, iMaxIter, iTimeStep, iDiffuCoeff, iField, ioField);
  else
    PreconditionedConjugateGradientSolve(iFieldID, iMaxIter, iTimeStep, iDiffuCoeff, iField, ioField);
}


// Solve linear system with Jacobi Preconditioned Conjugate Gradient approach
// Page 51 https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// https://github.com/awesson/stable-fluids/blob/master/linearSolver.cpp
// https://en.wikipedia.org/wiki/Conjugate_gradient_method
void CompuFluidDyna::PreconditionedConjugateGradientSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                                                          const std::vector<std::vector<std::vector<float>>>& iField,
                                                          std::vector<std::vector<std::vector<float>>>& ioField) {
  // Allocate and initialize variables and fields
  int stagnationCount= 0;
  float normSqrOld= 0.0f, normSqrChg= 0.0f;
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> dField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> sField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // r= b - A x
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  // d= M^-1 r
  ImplicitFieldLaplacianDiagPrecond(iFieldID, iTimeStep, iDiffuCoeff, rField, dField);
  // errNew= r^T d
  const float errBeg= ImplicitFieldDotProd(rField, dField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[PlotSolve_______].B()) {
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
    ApplyBC(iFieldID, t0Field);
    ImplicitFieldSub(iField, t0Field, t1Field);
    const float errTmp= ImplicitFieldDotProd(t1Field, t1Field);
    D.Plot[D.Plot.size() - 1].val.push_back(errTmp);
  }
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= 0.0f) break;
    if (errNew / errBeg <= D.UI[SolvTolRelResid_].F()) break;
    if (normSqrOld > 0.0f && normSqrChg / normSqrOld <= D.UI[SolvTolChgPrev__].F()) stagnationCount++;
    else stagnationCount= 0;
    if (stagnationCount > std::max(D.UI[SolvTolStagna___].I(), 0)) break;
    // q= A d
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, dField, qField);
    // alpha= errNew / (d^T q)
    const float denom= ImplicitFieldDotProd(dField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // temp= x + alpha d
    ImplicitFieldScale(alpha, dField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ApplyBC(iFieldID, t1Field);
    // Compute relative solution change
    if (D.UI[SolvTolChgPrev__].F() > 0.0f) {
      normSqrChg= ImplicitFieldDotProd(t0Field, t0Field);
      normSqrOld= ImplicitFieldDotProd(ioField, ioField);
    }
    // x= temp
    ioField= t1Field;
    // Compute residual
    if ((idxIter + 1) % std::max(D.UI[SolvDriftReset__].I(), 1) == 0) {
      // r= b - A x
      ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
      ApplyBC(iFieldID, t0Field);
      ImplicitFieldSub(iField, t0Field, rField);
    }
    else {
      // r= r - alpha q
      ImplicitFieldScale(alpha, qField, t0Field);
      ImplicitFieldSub(rField, t0Field, t1Field);
      rField= t1Field;
    }
    // s= M^-1 r
    ImplicitFieldLaplacianDiagPrecond(iFieldID, iTimeStep, iDiffuCoeff, rField, sField);
    // errNew= r^T s
    const float errOld= errNew;
    errNew= ImplicitFieldDotProd(rField, sField);
    // Error plot
    if (D.UI[PlotSolve_______].B()) {
      ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
      ApplyBC(iFieldID, t0Field);
      ImplicitFieldSub(iField, t0Field, t1Field);
      const float errTmp= ImplicitFieldDotProd(t1Field, t1Field);
      D.Plot[D.Plot.size() - 1].val.push_back(errTmp);
    }
    // d= s + (errNew / errOld) * d
    ImplicitFieldScale(errNew / errOld, dField, t0Field);
    ImplicitFieldAdd(sField, t0Field, dField);
  }
  // Save residual field for monitoring
  if (D.UI[PlotSolve_______].B() && D.UI[PlotSolve_______].I() == (int)D.Plot.size()) Dum0= rField;
}


// Solve linear system with Conjugate Gradient approach
// Page 50 https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// https://github.com/awesson/stable-fluids/blob/master/linearSolver.cpp
// https://en.wikipedia.org/wiki/Conjugate_gradient_method
void CompuFluidDyna::ConjugateGradientSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                                            const std::vector<std::vector<std::vector<float>>>& iField,
                                            std::vector<std::vector<std::vector<float>>>& ioField) {
  // Allocate and initialize variables and fields
  int stagnationCount= 0;
  float normSqrOld= 0.0f, normSqrChg= 0.0f;
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> dField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // r= b - A x
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  // d= r
  dField= rField;
  // errNew= r^T r
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[PlotSolve_______].B()) D.Plot[D.Plot.size() - 1].val.push_back(errNew);
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= 0.0f) break;
    if (errNew / errBeg <= D.UI[SolvTolRelResid_].F()) break;
    if (normSqrOld > 0.0f && normSqrChg / normSqrOld <= D.UI[SolvTolChgPrev__].F()) stagnationCount++;
    else stagnationCount= 0;
    if (stagnationCount > std::max(D.UI[SolvTolStagna___].I(), 0)) break;
    // q= A d
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, dField, qField);
    // alpha= errNew / (d^T q)
    const float denom= ImplicitFieldDotProd(dField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // temp= x + alpha d
    ImplicitFieldScale(alpha, dField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ApplyBC(iFieldID, t1Field);
    // Compute relative solution change
    if (D.UI[SolvTolChgPrev__].F() > 0.0f) {
      normSqrChg= ImplicitFieldDotProd(t0Field, t0Field);
      normSqrOld= ImplicitFieldDotProd(ioField, ioField);
    }
    // x= temp
    ioField= t1Field;
    // Compute residual
    if ((idxIter + 1) % std::max(D.UI[SolvDriftReset__].I(), 1) == 0) {
      // r= b - A x
      ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
      ApplyBC(iFieldID, t0Field);
      ImplicitFieldSub(iField, t0Field, rField);
    }
    else {
      // r= r - alpha q
      ImplicitFieldScale(alpha, qField, t0Field);
      ImplicitFieldSub(rField, t0Field, t1Field);
      rField= t1Field;
    }
    // errNew= r^T r
    const float errOld= errNew;
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[PlotSolve_______].B()) D.Plot[D.Plot.size() - 1].val.push_back(errNew);
    // d= r + (errNew / errOld) * d
    ImplicitFieldScale(errNew / errOld, dField, t0Field);
    ImplicitFieldAdd(rField, t0Field, dField);
  }
  // Save residual field for monitoring
  if (D.UI[PlotSolve_______].B() && D.UI[PlotSolve_______].I() == (int)D.Plot.size()) Dum0= rField;
}


// Solve linear system with Gradient descent approach
// Page 49 https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
void CompuFluidDyna::GradientDescentSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                                          const std::vector<std::vector<std::vector<float>>>& iField,
                                          std::vector<std::vector<std::vector<float>>>& ioField) {
  // Allocate and initialize variables and fields
  int stagnationCount= 0;
  float normSqrOld= 0.0f, normSqrChg= 0.0f;
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // r= b - A x
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  // errNew= r^T r
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[PlotSolve_______].B()) D.Plot[D.Plot.size() - 1].val.push_back(errNew);
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= 0.0f) break;
    if (errNew / errBeg <= D.UI[SolvTolRelResid_].F()) break;
    if (normSqrOld > 0.0f && normSqrChg / normSqrOld <= D.UI[SolvTolChgPrev__].F()) stagnationCount++;
    else stagnationCount= 0;
    if (stagnationCount > std::max(D.UI[SolvTolStagna___].I(), 0)) break;
    // q= A r
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, rField, qField);
    // alpha= errNew / (r^T q)
    const float denom= ImplicitFieldDotProd(rField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // temp= x + alpha r
    ImplicitFieldScale(alpha, rField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ApplyBC(iFieldID, t1Field);
    // Compute relative solution change
    if (D.UI[SolvTolChgPrev__].F() > 0.0f) {
      normSqrChg= ImplicitFieldDotProd(t0Field, t0Field);
      normSqrOld= ImplicitFieldDotProd(ioField, ioField);
    }
    // x= temp
    ioField= t1Field;
    // Compute residual
    if ((idxIter + 1) % std::max(D.UI[SolvDriftReset__].I(), 1) == 0) {
      // r= b - A x
      ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
      ApplyBC(iFieldID, t0Field);
      ImplicitFieldSub(iField, t0Field, rField);
    }
    else {
      // r= r - alpha q
      ImplicitFieldScale(alpha, qField, t0Field);
      ImplicitFieldSub(rField, t0Field, t1Field);
      rField= t1Field;
    }
    // errNew= r^T r
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[PlotSolve_______].B()) D.Plot[D.Plot.size() - 1].val.push_back(errNew);
  }
  // Save residual field for monitoring
  if (D.UI[PlotSolve_______].B() && D.UI[PlotSolve_______].I() == (int)D.Plot.size()) Dum0= rField;
}


// Solve linear system with an iterative Gauss Seidel scheme
// Solving each equation sequentially by cascading latest solution to next equation
// Run a forward and backward pass in parallel to avoid element ordering bias
// Apply successive overrelaxation coefficient to accelerate convergence
void CompuFluidDyna::GaussSeidelSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                                      const std::vector<std::vector<std::vector<float>>>& iField,
                                      std::vector<std::vector<std::vector<float>>>& ioField) {
  // TODO Fix Gauss Seidel behavior with porosity field
  // Allocate and initialize variables and fields
  int stagnationCount= 0;
  float normSqrOld= 0.0f, normSqrChg= 0.0f;
  const float voxSizeSqr= voxSize * voxSize;
  const float diffuVal= iDiffuCoeff * iTimeStep / voxSizeSqr;
  const float coeffOverrelax= 1.8f;
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<std::vector<float>>>> FieldT(2);
  // r= b - A x
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  // errNew= r^T r
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[PlotSolve_______].B()) D.Plot[D.Plot.size() - 1].val.push_back(errNew);
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= 0.0f) break;
    if (errNew / errBeg <= D.UI[SolvTolRelResid_].F()) break;
    if (normSqrOld > 0.0f && normSqrChg / normSqrOld <= D.UI[SolvTolChgPrev__].F()) stagnationCount++;
    else stagnationCount= 0;
    if (stagnationCount > std::max(D.UI[SolvTolStagna___].I(), 0)) break;
    // Initialize fields for forward and backward passes
    FieldT[0]= ioField;
    FieldT[1]= ioField;
    // Execute the two passes in parallel
#pragma omp parallel for if (D.UI[Multithread_____].B())
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
            if (D.UI[TestPoroLapla___].B()) {
              if (x - 1 >= 0) sum+= (1.0f - Poros[x - 1][y][z]) * xBCVal + Poros[x - 1][y][z] * FieldT[k][x - 1][y][z];
              if (x + 1 < nX) sum+= (1.0f - Poros[x + 1][y][z]) * xBCVal + Poros[x + 1][y][z] * FieldT[k][x + 1][y][z];
              if (y - 1 >= 0) sum+= (1.0f - Poros[x][y - 1][z]) * yBCVal + Poros[x][y - 1][z] * FieldT[k][x][y - 1][z];
              if (y + 1 < nY) sum+= (1.0f - Poros[x][y + 1][z]) * yBCVal + Poros[x][y + 1][z] * FieldT[k][x][y + 1][z];
              if (z - 1 >= 0) sum+= (1.0f - Poros[x][y][z - 1]) * zBCVal + Poros[x][y][z - 1] * FieldT[k][x][y][z - 1];
              if (z + 1 < nZ) sum+= (1.0f - Poros[x][y][z + 1]) * zBCVal + Poros[x][y][z + 1] * FieldT[k][x][y][z + 1];
            }
            else {
              if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : FieldT[k][x - 1][y][z];
              if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : FieldT[k][x + 1][y][z];
              if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : FieldT[k][x][y - 1][z];
              if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : FieldT[k][x][y + 1][z];
              if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : FieldT[k][x][y][z - 1];
              if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : FieldT[k][x][y][z + 1];
            }
            // Set new value according to coefficients and flags
            if (count > 0) {
              const float prevVal= FieldT[k][x][y][z];
              if (iFieldID == FieldID::IDPres)
                FieldT[k][x][y][z]= (voxSizeSqr * iField[x][y][z] + sum) / (float)count;
              else
                FieldT[k][x][y][z]= (iField[x][y][z] + diffuVal * sum) / (1.0f + (float)count * diffuVal);
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
          t1Field[x][y][z]= (FieldT[0][x][y][z] + FieldT[1][x][y][z]) / 2.0f;
    if (D.UI[SolvTolChgPrev__].F() > 0.0f) {
      ImplicitFieldSub(ioField, t1Field, t0Field);
      normSqrChg= ImplicitFieldDotProd(t0Field, t0Field);
      normSqrOld= ImplicitFieldDotProd(ioField, ioField);
    }
    ioField= t1Field;
    // r= b - A x
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuCoeff, ioField, t0Field);
    ApplyBC(iFieldID, t0Field);
    ImplicitFieldSub(iField, t0Field, rField);
    // errNew= r^T r
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[PlotSolve_______].B()) D.Plot[D.Plot.size() - 1].val.push_back(errNew);
  }
  // Save residual field for monitoring
  if (D.UI[PlotSolve_______].B() && D.UI[PlotSolve_______].I() == (int)D.Plot.size()) Dum0= rField;
}
