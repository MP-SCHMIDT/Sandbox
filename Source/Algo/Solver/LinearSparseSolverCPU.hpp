#pragma once

// Standard lib
#include <vector>


class LinearSparseSolverCPU
{
  private:
  bool isSetup;               // Flag to mark if the matrix was set up
  int nbDof;                  // Number of dof
  int nbDofCompact;           // Number of dof when skipped dofs are excluded
  std::vector<int> expandor;  // Array to map between dof indices with/without skipping
  std::vector<int> crsRow;    // Sparse matrix in CRS format
  std::vector<int> crsCol;    // Sparse matrix in CRS format
  std::vector<float> crsVal;  // Sparse matrix in CRS format
  std::vector<float> precond; // Diagonal preconditionner coefficients
  bool usePrecond;


  public:
  LinearSparseSolverCPU();

  // Sets the compacted number of free dofs and the CRS and expandor arrays
  // Uses the matrix triplets in ascending order (row major) and a flag for skipped DOFs
  void SetupMatrix(const std::vector<int>& iTripletRow,
                   const std::vector<int>& iTripletCol,
                   const std::vector<float>& iTripletVal,
                   const std::vector<bool>& iSkipDOF,
                   const bool iUsePrecond,
                   const int iVerboseLevel);

  // Solves for the unknown vector given the current CRS matrix, the RHS and an initial guess
  // The residual history is returned
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
  // Run the iterative conjugate gradient solver
  void RunIterativePCG(const int iMaxIter,
                       const float iTolResidual,
                       const std::vector<float>& iRHSCompact,
                       std::vector<float>& ioSolCompact,
                       std::vector<float>& oErrorHistory);
  void RunIterativeCG(const int iMaxIter,
                      const float iTolResidual,
                      const std::vector<float>& iRHSCompact,
                      std::vector<float>& ioSolCompact,
                      std::vector<float>& oErrorHistory);
};