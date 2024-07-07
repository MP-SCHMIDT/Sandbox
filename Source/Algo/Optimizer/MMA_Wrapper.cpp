
#include "MMA_Wrapper.hpp"

// Eigen lib
#include <Eigen/Core>
#include <Eigen/Dense>

// Algo headers
#include "Optimizer/MMA.hpp"


void MMA_Wrapper::mmasub(
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
    std::vector<double>& ioVar) {
  // Allocate and set the Eigen VectorXd
  Eigen::VectorXd xval= Eigen::VectorXd::Map(ioVar.data(), iNbVar);
  Eigen::VectorXd xmin= Eigen::VectorXd::Map(iVarMin.data(), iNbVar);
  Eigen::VectorXd xmax= Eigen::VectorXd::Map(iVarMax.data(), iNbVar);
  Eigen::VectorXd df0dx= Eigen::VectorXd::Map(iObjGrad.data(), iNbVar);
  Eigen::VectorXd fval= Eigen::VectorXd::Map(iConstrVal.data(), iNbConstr);
  Eigen::VectorXd xold1= Eigen::VectorXd::Map(ioVarOld.data(), iNbVar);
  Eigen::VectorXd xold2= Eigen::VectorXd::Map(ioVarOldOld.data(), iNbVar);
  Eigen::VectorXd low= Eigen::VectorXd::Map(ioAsympLow.data(), iNbVar);
  Eigen::VectorXd upp= Eigen::VectorXd::Map(ioAsympUpp.data(), iNbVar);

  // Allocate and set the Eigen MatrixXd
  Eigen::MatrixXd dfdx(iNbConstr, iNbVar);
  for (int idxObj= 0; idxObj < iNbConstr; idxObj++)
    for (int idxCAD= 0; idxCAD < iNbVar; idxCAD++)
      dfdx(idxObj, idxCAD)= iConstrGrad[idxObj][idxCAD];

  // Call the MMA optimizer
  Eigen::VectorXd ox(iNbVar);
  MMA::mmasub(iNbConstr, iNbVar, iIdxIter, xval, xmin, xmax, df0dx, fval, dfdx, xold1, xold2, low, upp, ox);

  // Retreive the updated data
  for (int idxCAD= 0; idxCAD < iNbVar; idxCAD++) {
    ioVarOld[idxCAD]= xold1(idxCAD);
    ioVarOldOld[idxCAD]= xold2(idxCAD);
    ioAsympLow[idxCAD]= low(idxCAD);
    ioAsympUpp[idxCAD]= upp(idxCAD);
    ioVar[idxCAD]= ox(idxCAD);
  }
}
