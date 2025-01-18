#include "LinearSparseSolverCPU.hpp"


// Standard lib
#include <array>
#include <cassert>
#include <cmath>
#include <vector>
// TODO include cassert and swap all throws throughout the code

// Algo headers
#include "Util/Timer.hpp"


LinearSparseSolverCPU::LinearSparseSolverCPU() {
  isSetup= false;
}


void LinearSparseSolverCPU::SetupMatrix(const std::vector<int>& iTripletRow,
                                        const std::vector<int>& iTripletCol,
                                        const std::vector<float>& iTripletVal,
                                        const std::vector<bool>& iSkipDOF,
                                        const bool iUsePrecond,
                                        const int iVerboseLevel) {
  if (iVerboseLevel >= 3) Timer::PushTimer();

  nbDof= (int)iSkipDOF.size();
  const int nbTriplets= (int)iTripletVal.size();

  // Check triplets size ordering and non-zeroness
  assert(iTripletRow.size() == iTripletVal.size());
  assert(iTripletCol.size() == iTripletVal.size());
  for (int k= 0; k < nbTriplets; k++) {
    if (k+1 < nbTriplets) assert(iTripletRow[k] < iTripletRow[k+1] || (iTripletRow[k] == iTripletRow[k+1] && iTripletCol[k] < iTripletCol[k+1]));
    assert(iTripletVal[k] != 0.0f);
  }

  // Build the compactor expandor vectors
  std::vector<int> compactor(nbDof, -1);
  expandor= std::vector<int>(nbDof, -1);
  nbDofCompact= 0;
  for (int k= 0; k < nbDof; k++) {
    if (iSkipDOF[k]) continue;
    compactor[k]= nbDofCompact;
    expandor[nbDofCompact]= k;
    nbDofCompact++;
  }

  // Compactify the nnz row and column indices to ignore inactive elements
  std::vector<int> cooRow(nbTriplets);
  crsCol= std::vector<int>(nbTriplets);
  crsVal= iTripletVal;
  for (int k= 0; k < nbTriplets; k++) {
    cooRow[k]= compactor[iTripletRow[k]];
    crsCol[k]= compactor[iTripletCol[k]];
  }

  // Compress the row indices to encode in Compressed Row Storage (CRS) format
  int idxCRSNNZ= 0;
  crsRow= std::vector<int>(nbDofCompact+1);
  for (int idxRow= 0; idxRow < nbDofCompact; idxRow++) {
    crsRow[idxRow]= idxCRSNNZ;
    while (idxCRSNNZ < nbTriplets && cooRow[idxCRSNNZ] == idxRow)
      idxCRSNNZ++;
  }
  crsRow[nbDofCompact]= idxCRSNNZ;

  // Save the inverse of diagonal coefficients to be used as preconditionner
  if (iUsePrecond) {
    usePrecond= true;
    precond= std::vector<float>(nbDofCompact, 1.0f);
    #pragma omp parallel for
    for (int k= 0; k < nbDofCompact; k++)
      for (int idxNNZ= crsRow[k]; idxNNZ < crsRow[k + 1]; idxNNZ++)
        if (k == crsCol[idxNNZ] && crsVal[idxNNZ] != 0.0f)
          precond[k]= 1.0f / crsVal[idxNNZ];
  }
  else {
    usePrecond= false;
    precond.clear();
  }

  isSetup= true;

  if (iVerboseLevel >= 3) {
    printf("%5.2f setupT ", Timer::PopTimer());
    fflush(stdout);
  }
}


std::vector<float> LinearSparseSolverCPU::SolveProblem(const std::vector<float>& iRHS,
                                                       std::vector<float>& ioSolution,
                                                       const int iMaxIter,
                                                       const float iTolResidual,
                                                       const int iVerboseLevel) {
  if (iVerboseLevel >= 3) Timer::PushTimer();

  // Check inputs
  assert(isSetup);
  assert(nbDof == (int)iRHS.size());
  assert(nbDof == (int)ioSolution.size());

  // Create compacted rhs and solution arrays
  std::vector<float> rhsCompact(nbDofCompact);
  std::vector<float> solCompact(nbDofCompact);
  for (int k= 0; k < nbDofCompact; k++) {
    rhsCompact[k]= iRHS[expandor[k]];
    solCompact[k]= ioSolution[expandor[k]];
  }

  // Run the solver
  std::vector<float> errorHistory;
  if (usePrecond) RunIterativePCG(iMaxIter, iTolResidual, rhsCompact, solCompact, errorHistory);
  else            RunIterativeCG(iMaxIter, iTolResidual, rhsCompact, solCompact, errorHistory);

  // Read compacted solution array
  for (int k= 0; k < nbDofCompact; k++)
    ioSolution[expandor[k]]= solCompact[k];

  // Print solve info
  if (iVerboseLevel >= 3) {
    printf("%d dofs %5d itSolv %3.2e resid %5.2f solvCPUT ", nbDofCompact, int(errorHistory.size())-1, errorHistory[errorHistory.size()-1], Timer::PopTimer());
    if (int(errorHistory.size())-1 >= iMaxIter) printf("[Diverged] ");
    fflush(stdout);
  }

  return errorHistory;
}


std::vector<float> LinearSparseSolverCPU::SolveProblem(const std::vector<double>& iRHS,
                                                       std::vector<double>& ioSolution,
                                                       const int iMaxIter,
                                                       const float iTolResidual,
                                                       const int iVerboseLevel) {
  if (iVerboseLevel >= 3) Timer::PushTimer();

  // Check inputs
  assert(isSetup);
  assert(nbDof == (int)iRHS.size());
  assert(nbDof == (int)ioSolution.size());

  // Create compacted rhs and solution arrays
  std::vector<float> rhsCompact(nbDofCompact);
  std::vector<float> solCompact(nbDofCompact);
  for (int k= 0; k < nbDofCompact; k++) {
    rhsCompact[k]= (float)iRHS[expandor[k]];
    solCompact[k]= (float)ioSolution[expandor[k]];
  }

  // Run the solver
  std::vector<float> errorHistory;
  if (usePrecond) RunIterativePCG(iMaxIter, iTolResidual, rhsCompact, solCompact, errorHistory);
  else            RunIterativeCG(iMaxIter, iTolResidual, rhsCompact, solCompact, errorHistory);

  // Read compacted solution array
  for (int k= 0; k < nbDofCompact; k++)
    ioSolution[expandor[k]]= (double)solCompact[k];

  // Print solve info
  if (iVerboseLevel >= 3) {
    printf("%d dofs %5d itSolv %3.2e resid %5.2f solvCPUT ", nbDofCompact, int(errorHistory.size())-1, errorHistory[errorHistory.size()-1], Timer::PopTimer());
    if (int(errorHistory.size())-1 >= iMaxIter) printf("[Diverged] ");
    fflush(stdout);
  }

  return errorHistory;
}


void LinearSparseSolverCPU::RunIterativePCG(const int iMaxIter,
                                            const float iTolResidual,
                                            const std::vector<float>& iRHSCompact,
                                            std::vector<float>& ioSolCompact,
                                            std::vector<float>& oErrorHistory) {
  oErrorHistory.clear();

  constexpr int sizePartial= 64;
  const int nbPartial= nbDofCompact / sizePartial + ((nbDofCompact % sizePartial == 0) ? 0 : 1);
  std::vector<float> sumPartial(nbPartial, 0.0f);

  // Compute RHS squared norm used by stopping criterion
  std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
  #pragma omp parallel for
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
      sumPartial[idxPartial]+= iRHSCompact[k] * iRHSCompact[k];
  float rhsNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    rhsNorm2+= sumPartial[idxPartial];

  // Check exit condition
  if (rhsNorm2 == 0.0f) {
    ioSolCompact= std::vector<float>(nbDofCompact, 0.0f);
    oErrorHistory.push_back(0.0f);
    return;
  }

  // Compute the residual vector    r= b - A x
  std::vector<float> rField(nbDofCompact, 0.0f);
  #pragma omp parallel for
  for (int k= 0; k < nbDofCompact; k++) {
    for (int idxNNZ= crsRow[k]; idxNNZ < crsRow[k + 1]; idxNNZ++)
      rField[k]+= crsVal[idxNNZ] * ioSolCompact[crsCol[idxNNZ]];
    rField[k]= iRHSCompact[k] - rField[k];
  }

  // Compute residual squared norm
  std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
  #pragma omp parallel for
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
        sumPartial[idxPartial]+= rField[k] * rField[k];
  float residNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    residNorm2+= sumPartial[idxPartial];

  // Check exit condition
  oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
  const float residNorm2Tol= iTolResidual * iTolResidual * rhsNorm2;
  if (residNorm2 <= residNorm2Tol) return;

  // Prepare additional fields    d= M^-1 r
  std::vector<float> qField(nbDofCompact, 0.0f);
  std::vector<float> dField(nbDofCompact, 0.0f);
  std::vector<float> sField(nbDofCompact, 0.0f);
  for (int k= 0; k < nbDofCompact; k++)
    dField[k]= precond[k] * rField[k];

  // Compute the error    errNew= r^T d
  std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
  #pragma omp parallel for
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
      sumPartial[idxPartial]+= rField[k] * dField[k];
  float errNew= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    errNew+= sumPartial[idxPartial];

  // Iterate to solve
  int idxIter;
  for (idxIter= 0; idxIter < iMaxIter; idxIter++) {

    // Multiply modified residual with matrix   q= A d
    #pragma omp parallel for
    for (int k= 0; k < nbDofCompact; k++) {
      qField[k]= 0.0f;
      for (int idxNNZ= crsRow[k]; idxNNZ < crsRow[k + 1]; idxNNZ++)
        qField[k]+= crsVal[idxNNZ] * dField[crsCol[idxNNZ]];
    }

    // Line search distance    alpha= errNew / (d^T q)
    std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
    #pragma omp parallel for
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
        sumPartial[idxPartial]+= dField[k] * qField[k];
    float denom= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      denom+= sumPartial[idxPartial];
    assert(denom != 0.0f);
    const float alpha= errNew / denom;

    // Update solution and residual
    // x ⇐ x + alpha d
    // r ⇐ r - alpha q
    for (int k= 0; k < nbDofCompact; k++) {
      ioSolCompact[k]+= alpha * dField[k];
      rField[k]-= alpha * qField[k];
    }

    // Update preconditioned residual
    // s ⇐ M^-1 r
    for (int k= 0; k < nbDofCompact; k++)
      sField[k]= precond[k] * rField[k];

    // Compute residual squared norm
    std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
    #pragma omp parallel for
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
        sumPartial[idxPartial]+= rField[k] * rField[k];
    residNorm2= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      residNorm2+= sumPartial[idxPartial];

    // Check exit condition
    oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
    if (residNorm2 <= residNorm2Tol) return;

    // Compute error    errNew= r^T s
    const float errOld= errNew;
    std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
    #pragma omp parallel for
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
        sumPartial[idxPartial]+= rField[k] * sField[k];
    errNew= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      errNew+= sumPartial[idxPartial];

    // Compute Gram-Schmidt value for new search direction
    // d= s + (errNew / errOld) * d
    assert(errOld != 0.0f);
    const float beta= errNew / errOld;
    for (int k= 0; k < nbDofCompact; k++)
      dField[k]= sField[k] + beta * dField[k];
  }
}


void LinearSparseSolverCPU::RunIterativeCG(const int iMaxIter,
                                           const float iTolResidual,
                                           const std::vector<float>& iRHSCompact,
                                           std::vector<float>& ioSolCompact,
                                           std::vector<float>& oErrorHistory) {
  oErrorHistory.clear();

  constexpr int sizePartial= 64;
  const int nbPartial= nbDofCompact / sizePartial + ((nbDofCompact % sizePartial == 0) ? 0 : 1);
  std::vector<float> sumPartial(nbPartial, 0.0f);

  // Compute RHS squared norm used by stopping criterion
  std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
  #pragma omp parallel for
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
      sumPartial[idxPartial]+= iRHSCompact[k] * iRHSCompact[k];
  float rhsNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    rhsNorm2+= sumPartial[idxPartial];

  // Check exit condition
  if (rhsNorm2 == 0.0f) {
    ioSolCompact= std::vector<float>(nbDofCompact, 0.0f);
    oErrorHistory.push_back(0.0f);
    return;
  }

  // Compute the residual vector    r= b - A x
  std::vector<float> rField(nbDofCompact, 0.0f);
  #pragma omp parallel for
  for (int k= 0; k < nbDofCompact; k++) {
    for (int idxNNZ= crsRow[k]; idxNNZ < crsRow[k + 1]; idxNNZ++)
      rField[k]+= crsVal[idxNNZ] * ioSolCompact[crsCol[idxNNZ]];
    rField[k]= iRHSCompact[k] - rField[k];
  }

  // Compute residual squared norm
  std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
  #pragma omp parallel for
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
      sumPartial[idxPartial]+= rField[k] * rField[k];
  float residNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    residNorm2+= sumPartial[idxPartial];

  // Check exit condition
  oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
  const float residNorm2Tol= iTolResidual * iTolResidual * rhsNorm2;
  if (residNorm2 <= residNorm2Tol) return;

  // Prepare additional fields    d= r
  std::vector<float> qField(nbDofCompact, 0.0f);
  std::vector<float> dField(nbDofCompact, 0.0f);
  dField= rField;

  // Compute the error    errNew= r^T r
  float errNew= residNorm2;

  // Iterate to solve
  int idxIter;
  for (idxIter= 0; idxIter < iMaxIter; idxIter++) {

    // Multiply modified residual with matrix   q= A d
    #pragma omp parallel for
    for (int k= 0; k < nbDofCompact; k++) {
      qField[k]= 0.0f;
      for (int idxNNZ= crsRow[k]; idxNNZ < crsRow[k + 1]; idxNNZ++)
        qField[k]+= crsVal[idxNNZ] * dField[crsCol[idxNNZ]];
    }

    // Line search distance    alpha= errNew / (d^T q)
    std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
    #pragma omp parallel for
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
        sumPartial[idxPartial]+= dField[k] * qField[k];
    float denom= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      denom+= sumPartial[idxPartial];
    assert(denom != 0.0f);
    const float alpha= errNew / denom;

    // Update solution and residual
    // x ⇐ x + alpha d
    // r ⇐ r - alpha q
    for (int k= 0; k < nbDofCompact; k++) {
      ioSolCompact[k]+= alpha * dField[k];
      rField[k]-= alpha * qField[k];
    }

    // Compute residual squared norm
    std::fill(sumPartial.begin(), sumPartial.end(), 0.0f);
    #pragma omp parallel for
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      for (int k= idxPartial*sizePartial; k < std::min((idxPartial+1)*sizePartial, nbDofCompact); k++)
        sumPartial[idxPartial]+= rField[k] * rField[k];
    residNorm2= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      residNorm2+= sumPartial[idxPartial];

    // Check exit condition
    oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
    if (residNorm2 <= residNorm2Tol) return;

    // Compute error    errNew= r^T r
    const float errOld= errNew;
    errNew= residNorm2;

    // Compute Gram-Schmidt value for new search direction
    // d= r + (errNew / errOld) * d
    assert(errOld != 0.0f);
    const float beta= errNew / errOld;
    for (int k= 0; k < nbDofCompact; k++)
      dField[k]= rField[k] + beta * dField[k];
  }
}
