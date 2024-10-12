#pragma once

// Standard lib
#include <vector>


class LinearSparseSolverEigen
{
  public:
  static void IterativeSolver(
      std::vector<int> const& iARow,
      std::vector<int> const& iACol,
      std::vector<float> const& iAVal,
      std::vector<double> const& iB,
      std::vector<bool> const& iDLock,
      std::vector<double>& ioX,
      double const iSolveResidual,
      int const iSolveMaxCgIter,
      int const iVerboseLevel);
};
