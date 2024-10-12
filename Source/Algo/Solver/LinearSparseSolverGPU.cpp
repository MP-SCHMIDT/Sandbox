#include "LinearSparseSolverGPU.hpp"


// Standard lib
#include <array>
#include <cassert>
#include <cmath>
#include <vector>

// OpenCL lib and wrapper (hidden by pImpl idiom to avoid polluting the rest of the codebase)
#include "OpenCL_Wrapper/opencl.hpp" // Note: also includes the wrapper's ugly utilities.hpp header

// Algo headers
#include "Util/Timer.hpp"
#include "Solver/LinearSparseSolverGPU_Kernel.hpp"


class LinearSparseSolverGPU::Impl { // Definition of the implementation class
  private:
  bool isSetup;                    // Flag to mark if the matrix was set up
  int nbDof;                       // Number of dof
  int nbDofCompact;                // Number of dof when skipped dofs are excluded
  std::vector<int> expandor;       // Array to map between dof indices with/without skipping
  bool usePrecond;                 // Flag to mark whether the diagonal preconditionner is ready and in use

  struct {
    Device device;
    Kernel kernelMatrixMult;
    Kernel kernelVecMult;
    Kernel kernelAddScalarMult;
    Kernel kernelArrayCopy;
    Kernel kernelDotProdPartialReduc;
    Memory<int> arrMatRowStart;
    Memory<int> arrMatCol;
    Memory<float> arrMatVal;
    Memory<float> arrPrecond;
    Memory<float> arrRHS;
    Memory<float> arrSol;
    Memory<float> arrResid;
    Memory<float> arrQ;
    Memory<float> arrD;
    Memory<float> arrS;
    Memory<float> arrPartial;
  } OCL;

  public:
  Impl() {
      isSetup= false;
  }

  void SetupMatrix(const std::vector<int>& iTripletRow,
                                          const std::vector<int>& iTripletCol,
                                          const std::vector<float>& iTripletVal,
                                          const std::vector<bool>& iSkipDOF,
                                          const bool iUsePrecond,
                                          const int iVerboseLevel);


  std::vector<float> SolveProblem(const std::vector<float>& iRHS,
                                  std::vector<float>& ioSolution,
                                  const int iMaxIter,
                                  const float iTolResidual,
                                  const int iVerboseLevel);
  std::vector<float> SolveProblem(const std::vector<double>& iRHS,
                                  std::vector<double>& ioSolution,
                                  const int iMaxIter,
                                  const float iTolResidual,
                                  const int iVerboseLevel);

  private:
  void RunIterativePCG(const int iMaxIter, const float iTolResidual, std::vector<float>& oErrorHistory);
  void RunIterativeCG(const int iMaxIter, const float iTolResidual, std::vector<float>& oErrorHistory);

};


// Constructor for the Header class
LinearSparseSolverGPU::LinearSparseSolverGPU() : pImpl(std::make_unique<Impl>()) {}

// Destructor for the Header class
LinearSparseSolverGPU::~LinearSparseSolverGPU() = default; // Default destructor for pImpl

// Forward call of impl function
void LinearSparseSolverGPU::SetupMatrix(const std::vector<int>& iTripletRow,
                                        const std::vector<int>& iTripletCol,
                                        const std::vector<float>& iTripletVal,
                                        const std::vector<bool>& iSkipDOF,
                                        const bool iUsePrecond,
                                        const int iVerboseLevel) {
  pImpl->SetupMatrix(iTripletRow, iTripletCol, iTripletVal, iSkipDOF, iUsePrecond, iVerboseLevel); // Delegate the run function to the implementation
}

// Forward call of impl function
std::vector<float> LinearSparseSolverGPU::SolveProblem(const std::vector<float>& iRHS,
                                                       std::vector<float>& ioSolution,
                                                       const int iMaxIter,
                                                       const float iTolResidual,
                                                       const int iVerboseLevel) {
  return pImpl->SolveProblem(iRHS, ioSolution, iMaxIter, iTolResidual, iVerboseLevel); // Delegate the run function to the implementation
}
std::vector<float> LinearSparseSolverGPU::SolveProblem(const std::vector<double>& iRHS,
                                                       std::vector<double>& ioSolution,
                                                       const int iMaxIter,
                                                       const float iTolResidual,
                                                       const int iVerboseLevel) {
  return pImpl->SolveProblem(iRHS, ioSolution, iMaxIter, iTolResidual, iVerboseLevel); // Delegate the run function to the implementation
}


void LinearSparseSolverGPU::Impl::SetupMatrix(const std::vector<int>& iTripletRow,
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

  // Create the computation device
  OCL.device= Device((select_device_with_most_flops(get_devices(false))), LinearSparseSolverGPU_Kernel::get_opencl_c_code(), false);

  // Allocate the memory shared by host and device
  OCL.arrMatRowStart= Memory<int>(OCL.device, nbDofCompact + 1u, 1u, true, true, 0);
  OCL.arrMatCol=      Memory<int>(OCL.device, nbTriplets, 1u, true, true, 0);
  OCL.arrMatVal=      Memory<float>(OCL.device, nbTriplets, 1u, true, true, 0.0f);
  OCL.arrPrecond=     Memory<float>(OCL.device, nbDofCompact, 1u, true, true, 1.0f);
  OCL.arrRHS=         Memory<float>(OCL.device, nbDofCompact, 1u, true, true, 0.0f);
  OCL.arrSol=         Memory<float>(OCL.device, nbDofCompact, 1u, true, true, 0.0f);
  OCL.arrResid=       Memory<float>(OCL.device, nbDofCompact, 1u, false, true, 0.0f);
  OCL.arrQ=           Memory<float>(OCL.device, nbDofCompact, 1u, false, true, 0.0f);
  OCL.arrD=           Memory<float>(OCL.device, nbDofCompact, 1u, false, true, 0.0f);
  OCL.arrS=           Memory<float>(OCL.device, nbDofCompact, 1u, false, true, 0.0f);
  OCL.arrPartial=     Memory<float>(OCL.device, nbDofCompact / WORKGROUP_SIZE + 1u, 1u, true, true, 0.0f);

  // Create the OpenCL kernels
  OCL.kernelMatrixMult=          Kernel(OCL.device, nbDofCompact, "kernel_MatrixMult", OCL.arrMatRowStart, OCL.arrMatCol, OCL.arrMatVal, OCL.arrSol, OCL.arrS, nbDofCompact);
  OCL.kernelVecMult=             Kernel(OCL.device, nbDofCompact, "kernel_VecMult", OCL.arrPrecond, OCL.arrResid, OCL.arrD, nbDofCompact);
  OCL.kernelAddScalarMult=       Kernel(OCL.device, nbDofCompact, "kernel_AddScalarMult", OCL.arrRHS, -1.0f, OCL.arrS, OCL.arrResid, nbDofCompact);
  OCL.kernelArrayCopy=           Kernel(OCL.device, nbDofCompact, "kernel_ArrayCopy", OCL.arrResid, OCL.arrD, nbDofCompact);
  OCL.kernelDotProdPartialReduc= Kernel(OCL.device, nbDofCompact, "kernel_DotProdPartialReduc", OCL.arrResid, OCL.arrResid, OCL.arrPartial, nbDofCompact);

  // Compactify the nnz row and column indices to ignore inactive elements
  std::vector<int> cooRow(nbTriplets);
  for (int k= 0; k < nbTriplets; k++) {
    cooRow[k]= compactor[iTripletRow[k]];
    OCL.arrMatCol[k]= compactor[iTripletCol[k]];
    OCL.arrMatVal[k]= iTripletVal[k];
  }

  // Compress the row indices to encode in Compressed Row Storage (CRS) format
  int idxCRSNNZ= 0;
  for (int idxRow= 0; idxRow < nbDofCompact; idxRow++) {
    OCL.arrMatRowStart[idxRow]= idxCRSNNZ;
    while (idxCRSNNZ < nbTriplets && cooRow[idxCRSNNZ] == idxRow)
      idxCRSNNZ++;
  }
  OCL.arrMatRowStart[nbDofCompact]= idxCRSNNZ;

  // Save the inverse of diagonal coefficients to be used as preconditionner
  if (iUsePrecond) {
    usePrecond= true;
    #pragma omp parallel for
    for (int k= 0; k < nbDofCompact; k++)
      for (int idxNNZ= OCL.arrMatRowStart[k]; idxNNZ < OCL.arrMatRowStart[k + 1]; idxNNZ++)
        if (k == OCL.arrMatCol[idxNNZ] && OCL.arrMatVal[idxNNZ] != 0.0f)
          OCL.arrPrecond[k]= 1.0f / OCL.arrMatVal[idxNNZ];
  }
  else {
    usePrecond= false;
  }

  // Send the sparse matrix data to the device
  OCL.arrMatRowStart.write_to_device();
  OCL.arrMatCol.write_to_device();
  OCL.arrMatVal.write_to_device();
  if (usePrecond) OCL.arrPrecond.write_to_device();

  isSetup= true;

  if (iVerboseLevel >= 3) {
    printf("%5.2f setupT ", Timer::PopTimer());
    fflush(stdout);
  }
}


std::vector<float> LinearSparseSolverGPU::Impl::SolveProblem(const std::vector<float>& iRHS,
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
  for (int k= 0; k < nbDofCompact; k++) {
    OCL.arrRHS[k]= iRHS[expandor[k]];
    OCL.arrSol[k]= ioSolution[expandor[k]];
  }
  OCL.arrRHS.write_to_device();
  OCL.arrSol.write_to_device();

  // Run the solver
  std::vector<float> errorHistory;
  if (usePrecond) RunIterativePCG(iMaxIter, iTolResidual, errorHistory);
  else            RunIterativeCG(iMaxIter, iTolResidual, errorHistory);

  // Read compacted solution array
  OCL.arrSol.read_from_device();
  for (int k= 0; k < nbDofCompact; k++)
    ioSolution[expandor[k]]= OCL.arrSol[k];

  // Print solve info
  if (iVerboseLevel >= 3) {
    printf("%d dofs %5d itSolv %3.2e resid %5.2f solvGPUT ", nbDofCompact, int(errorHistory.size())-1, errorHistory[errorHistory.size()-1], Timer::PopTimer());
    if (int(errorHistory.size())-1 >= iMaxIter) printf("[Diverged] ");
    fflush(stdout);
  }

  return errorHistory;
}


std::vector<float> LinearSparseSolverGPU::Impl::SolveProblem(const std::vector<double>& iRHS,
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
  for (int k= 0; k < nbDofCompact; k++) {
    OCL.arrRHS[k]= (float)iRHS[expandor[k]];
    OCL.arrSol[k]= (float)ioSolution[expandor[k]];
  }
  OCL.arrRHS.write_to_device();
  OCL.arrSol.write_to_device();

  // Run the solver
  std::vector<float> errorHistory;
  if (usePrecond) RunIterativePCG(iMaxIter, iTolResidual, errorHistory);
  else            RunIterativeCG(iMaxIter, iTolResidual, errorHistory);

  // Read compacted solution array
  OCL.arrSol.read_from_device();
  for (int k= 0; k < nbDofCompact; k++)
    ioSolution[expandor[k]]= (double)OCL.arrSol[k];

  // Print solve info
  if (iVerboseLevel >= 3) {
    printf("%d dofs %5d itSolv %3.2e resid %5.2f solvGPUT ", nbDofCompact, int(errorHistory.size())-1, errorHistory[errorHistory.size()-1], Timer::PopTimer());
    if (int(errorHistory.size())-1 >= iMaxIter) printf("[Diverged] ");
    fflush(stdout);
  }

  return errorHistory;
}


void LinearSparseSolverGPU::Impl::RunIterativePCG(const int iMaxIter,
                                                  const float iTolResidual,
                                                  std::vector<float>& oErrorHistory) {
  oErrorHistory.clear();
  const int nbPartial= (int)OCL.arrPartial.length();

  // Compute RHS squared norm used by stopping criterion
  OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrRHS).set_parameters(1u, OCL.arrRHS).run();
  OCL.arrPartial.read_from_device();
  float rhsNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    rhsNorm2+= OCL.arrPartial[idxPartial];

  // Check exit condition
  if (rhsNorm2 == 0.0f) {
    for (int k= 0; k < nbDofCompact; k++)
      OCL.arrSol[k]= 0.0f;
    OCL.arrSol.write_to_device();
    oErrorHistory.push_back(0.0f);
    return;
  }

  // Compute the residual vector    r= b - A x
  OCL.kernelMatrixMult.set_parameters(3u, OCL.arrSol).set_parameters(4u, OCL.arrS).run();
  OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrRHS).set_parameters(1u, -1.0f).set_parameters(2u, OCL.arrS).set_parameters(3u, OCL.arrResid).run();

  // Compute residual squared norm
  OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrResid).run();
  OCL.arrPartial.read_from_device();
  float residNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    residNorm2+= OCL.arrPartial[idxPartial];

  // Check exit condition
  oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
  const float residNorm2Tol= iTolResidual * iTolResidual * rhsNorm2;
  if (residNorm2 < residNorm2Tol) return;

  // Prepare additional fields    d= M^-1 r
  OCL.kernelVecMult.set_parameters(1u, OCL.arrResid).set_parameters(2u, OCL.arrD).run();

  // Compute the error    errNew= r^T d
  OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrD).run();
  OCL.arrPartial.read_from_device();
  float errNew= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    errNew+= OCL.arrPartial[idxPartial];

  // Iterate to solve
  int idxIter;
  for (idxIter= 0; idxIter < iMaxIter; idxIter++) {

    // Multiply modified residual with matrix   q= A d
    OCL.kernelMatrixMult.set_parameters(3u, OCL.arrD).set_parameters(4u, OCL.arrQ).run();

    // Line search distance    alpha= errNew / (d^T q)
    OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrD).set_parameters(1u, OCL.arrQ).run();
    OCL.arrPartial.read_from_device();
    float denom= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      denom+= OCL.arrPartial[idxPartial];
    assert(denom != 0.0f);
    const float alpha= errNew / denom;

    // Update solution and residual
    // x ⇐ x + alpha d
    // r ⇐ r - alpha q
    // s ⇐ M^-1 r
    OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrSol).set_parameters(1u, alpha).set_parameters(2u, OCL.arrD).set_parameters(3u, OCL.arrSol).run();
    OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrResid).set_parameters(1u, -alpha).set_parameters(2u, OCL.arrQ).set_parameters(3u, OCL.arrResid).run();

    // Update preconditioned residual
    // s ⇐ M^-1 r
    OCL.kernelVecMult.set_parameters(1u, OCL.arrResid).set_parameters(2u, OCL.arrS).run();

    // Compute residual squared norm
    OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrResid).run();
    OCL.arrPartial.read_from_device();
    residNorm2= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      residNorm2+= OCL.arrPartial[idxPartial];

    // Check exit condition
    oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
    if (residNorm2 < residNorm2Tol) return;

    // Compute error    errNew= r^T s
    const float errOld= errNew;
    OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrS).run();
    OCL.arrPartial.read_from_device();
    errNew= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      errNew+= OCL.arrPartial[idxPartial];

    // Compute Gram-Schmidt value for new search direction
    // d= s + (errNew / errOld) * d
    assert(errOld != 0.0f);
    const float beta= errNew / errOld;
    OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrS).set_parameters(1u, beta).set_parameters(2u, OCL.arrD).set_parameters(3u, OCL.arrD).run();
  }
}



void LinearSparseSolverGPU::Impl::RunIterativeCG(const int iMaxIter,
                                                 const float iTolResidual,
                                                 std::vector<float>& oErrorHistory) {
  oErrorHistory.clear();
  const int nbPartial= (int)OCL.arrPartial.length();

  // Compute RHS squared norm used by stopping criterion
  OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrRHS).set_parameters(1u, OCL.arrRHS).run();
  OCL.arrPartial.read_from_device();
  float rhsNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    rhsNorm2+= OCL.arrPartial[idxPartial];

  // Check exit condition
  if (rhsNorm2 == 0.0f) {
    for (int k= 0; k < nbDofCompact; k++)
      OCL.arrSol[k]= 0.0f;
    OCL.arrSol.write_to_device();
    oErrorHistory.push_back(0.0f);
    return;
  }

  // Compute the residual vector    r= b - A x
  OCL.kernelMatrixMult.set_parameters(3u, OCL.arrSol).set_parameters(4u, OCL.arrS).run();
  OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrRHS).set_parameters(1u, -1.0f).set_parameters(2u, OCL.arrS).set_parameters(3u, OCL.arrResid).run();

  // Compute residual squared norm
  OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrResid).run();
  OCL.arrPartial.read_from_device();
  float residNorm2= 0.0f;
  for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
    residNorm2+= OCL.arrPartial[idxPartial];

  // Check exit condition
  oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
  const float residNorm2Tol= iTolResidual * iTolResidual * rhsNorm2;
  if (residNorm2 < residNorm2Tol) return;

  // Prepare additional fields    d= r
  OCL.kernelArrayCopy.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrD).run();

  // Compute the error    errNew= r^T r
  float errNew= residNorm2;

  // Iterate to solve
  int idxIter;
  for (idxIter= 0; idxIter < iMaxIter; idxIter++) {

    // Multiply modified residual with matrix   q= A d
    OCL.kernelMatrixMult.set_parameters(3u, OCL.arrD).set_parameters(4u, OCL.arrQ).run();

    // Line search distance    alpha= errNew / (d^T q)
    OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrD).set_parameters(1u, OCL.arrQ).run();
    OCL.arrPartial.read_from_device();
    float denom= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      denom+= OCL.arrPartial[idxPartial];
    assert(denom != 0.0f);
    const float alpha= errNew / denom;

    // Update solution and residual
    // x ⇐ x + alpha d
    // r ⇐ r - alpha q
    OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrSol).set_parameters(1u, alpha).set_parameters(2u, OCL.arrD).set_parameters(3u, OCL.arrSol).run();
    OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrResid).set_parameters(1u, -alpha).set_parameters(2u, OCL.arrQ).set_parameters(3u, OCL.arrResid).run();

    // Compute residual squared norm
    OCL.kernelDotProdPartialReduc.set_parameters(0u, OCL.arrResid).set_parameters(1u, OCL.arrResid).run();
    OCL.arrPartial.read_from_device();
    residNorm2= 0.0f;
    for (int idxPartial= 0; idxPartial < nbPartial; idxPartial++)
      residNorm2+= OCL.arrPartial[idxPartial];

    // Check exit condition
    oErrorHistory.push_back(std::sqrt(residNorm2 / rhsNorm2));
    if (residNorm2 < residNorm2Tol) return;

    // Compute error    errNew= r^T r
    const float errOld= errNew;
    errNew= residNorm2;

    // Compute Gram-Schmidt value for new search direction
    // d= r + (errNew / errOld) * d
    assert(errOld != 0.0f);
    const float beta= errNew / errOld;
    OCL.kernelAddScalarMult.set_parameters(0u, OCL.arrResid).set_parameters(1u, beta).set_parameters(2u, OCL.arrD).set_parameters(3u, OCL.arrD).run();
  }
}
