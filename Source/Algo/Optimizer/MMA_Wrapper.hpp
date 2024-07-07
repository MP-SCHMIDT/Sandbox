#pragma once

// Standard lib
#include <vector>


class MMA_Wrapper
{
  public:
  // Wapper method for the MMA Optimizer in Source\Algo\Optimizer\MMA.hpp
  static void mmasub(
      int const iNbVar,
      int const iNbConstr,
      int const iIdxIter,
      std::vector<double> const& iVarMin,
      std::vector<double> const& iVarMax,
      std::vector<double> const& iObjGrad,
      std::vector<double> const& iConstrVal,
      std::vector<std::vector<double>> const& iConstrGrad,
      std::vector<double>& ioAsympLow,
      std::vector<double>& ioAsympUpp,
      std::vector<double>& ioVarOldOld,
      std::vector<double>& ioVarOld,
      std::vector<double>& ioVar);
};
