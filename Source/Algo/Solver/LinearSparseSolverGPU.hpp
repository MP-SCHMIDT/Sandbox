#pragma once

// Standard lib
#include <vector>
#include <memory> // Include for std::unique_ptr


class LinearSparseSolverGPU
{
  public:
  std::vector<float> errorHistory; // Residual error values of the last iterative solve

  LinearSparseSolverGPU();
  ~LinearSparseSolverGPU();

  // Sets the compacted number of free dofs and the CRS and expandor arrays
  // Uses the matrix triplets in ascending order (row major) and a flag for skipped DOFs
  void SetupMatrix(const std::vector<int>& iTripletRow,
                   const std::vector<int>& iTripletCol,
                   const std::vector<float>& iTripletVal,
                   const std::vector<bool>& iSkipDOF,
                   const bool iUsePrecond,
                   const bool iVerbose);

  // Solves for the unknown vector given the current CRS matrix, the RHS and an initial guess
  // The residual history is returned
  void SolveProblem(const std::vector<float>& iRHS,
                    std::vector<float>& ioSolution,
                    const int iMaxIter,
                    const float iTolResidual,
                    const bool iVerbose);
  void SolveProblem(const std::vector<double>& iRHS,
                    std::vector<double>& ioSolution,
                    const int iMaxIter,
                    const float iTolResidual,
                    const bool iVerbose);

  private:
  class Impl;                        // Forward declaration of the private implementation class
  std::unique_ptr<Impl> pImpl;       // Pointer to the implementation
};
