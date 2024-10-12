
#include "LinearSparseSolverEigen.hpp"


// Standard lib
#include <vector>

// Eigen lib
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

// Algo headers
#include "Util/Timer.hpp"


// #define TESTING_AND_DEBUG_MODE
#ifdef TESTING_AND_DEBUG_MODE
#pragma optimize("", off)
#endif
void LinearSparseSolverEigen::IterativeSolver(
    std::vector<int> const& iARow,
    std::vector<int> const& iACol,
    std::vector<float> const& iAVal,
    std::vector<double> const& iB,
    std::vector<bool> const& iDLock,
    std::vector<double>& ioX,
    double const iSolveResidual,
    int const iSolveMaxCgIter,
    int const iVerboseLevel) {
  if (iVerboseLevel >= 3) Timer::PushTimer();

  // Build the compactor expandor vectors
  int const nbDof= (int)iB.size();
  std::vector<int> compactor(nbDof, -1);
  std::vector<int> expandor(nbDof, -1);
  int nbDofCompact= 0;
  for (int k= 0; k < nbDof; k++) {
    if (iDLock[k]) continue;
    compactor[k]= nbDofCompact;
    expandor[nbDofCompact]= k;
    nbDofCompact++;
  }

  // Build the sparse matrix
  Eigen::SparseMatrix<float> A(nbDofCompact, nbDofCompact);
  {
    std::vector<Eigen::Triplet<float, int>> tripletList(iAVal.size());
    for (int k= 0; k < (int)iAVal.size(); k++)
      tripletList[k]= Eigen::Triplet<float, int>(compactor[iARow[k]], compactor[iACol[k]], iAVal[k]);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A.makeCompressed();
  }

  // Build the equation vectors
  Eigen::VectorXf XGuess(nbDofCompact);
  Eigen::VectorXf B(nbDofCompact);
  for (int k= 0; k < nbDofCompact; k++) {
    B[k]= float(iB[expandor[k]]);
    XGuess[k]= float(ioX[expandor[k]]);
  }

  // Setup the solver
  Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper> cg;
  cg.setTolerance(float(iSolveResidual));
  cg.setMaxIterations(Eigen::Index(iSolveMaxCgIter));
  cg.compute(A);

  if (iVerboseLevel >= 3) {
    printf("%5.2f setupT ", Timer::PopTimer());
    fflush(stdout);
  }

  // Solve the system
  if (iVerboseLevel >= 3) Timer::PushTimer();
  Eigen::VectorXf X= cg.solveWithGuess(B, XGuess);

  // Copy the result in the output vector
  ioX= std::vector<double>(nbDof, 0.0);
  for (int k= 0; k < nbDofCompact; k++) {
    ioX[expandor[k]]= double(X[k]);
  }
  if (iVerboseLevel >= 3) {
    printf("%d dofs %5d itSolv %3.2e resid %5.2f solvEigT ", nbDofCompact, (int)cg.iterations(), (double)cg.error(), Timer::PopTimer());
    fflush(stdout);
  }
}
#ifdef TESTING_AND_DEBUG_MODE
#pragma optimize("", on)
#endif
