#include "NonLinMMABench.hpp"


// Standard lib

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Eigen lib
#include <Eigen/Core>
#include <Eigen/Dense>

// Algo headers
#include "Optimizer/MMA.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
NonLinMMABench::NonLinMMABench() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void NonLinMMABench::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ResetOptim______", 0));
    D.UI.push_back(ParamUI("ForceSolution___", 0));
    D.UI.push_back(ParamUI("InitVal_________", 5.0));
    D.UI.push_back(ParamUI("MinVal__________", 1.e-6));
    D.UI.push_back(ParamUI("MaxVal__________", 1.e6));
    D.UI.push_back(ParamUI("MoveLimit_______", 0.5));
    D.UI.push_back(ParamUI("ConstraintVal___", 1.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool NonLinMMABench::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool NonLinMMABench::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void NonLinMMABench::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Allocate()\n");
}


// Refresh the project
void NonLinMMABench::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Refresh()\n");
}


// Handle UI parameter change
void NonLinMMABench::ParamChange() {
}


// Handle keypress
void NonLinMMABench::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  if (D.UI[VerboseLevel____].I() >= 5) printf("KeyPress()\n");

  // Allocate stuff
  static int itera= 0;
  const int nbConstraints= 1;
  const int nbElemOptimized= 5;

  static Eigen::VectorXd ix(nbElemOptimized);
  static Eigen::VectorXd ixmin(nbElemOptimized);
  static Eigen::VectorXd ixmax(nbElemOptimized);
  static Eigen::VectorXd df0dx(nbElemOptimized);
  static Eigen::VectorXd fval(nbConstraints);
  static Eigen::MatrixXd dfdx(nbConstraints, nbElemOptimized);
  static Eigen::VectorXd ox(nbElemOptimized);

  static Eigen::VectorXd xold1(nbElemOptimized);
  static Eigen::VectorXd xold2(nbElemOptimized);
  static Eigen::VectorXd low(nbElemOptimized);
  static Eigen::VectorXd upp(nbElemOptimized);

  // Force reset based on user request
  if (D.UI[ResetOptim______].I() > 0) {
    itera= 0;
  }

  // Initialize arrays
  if (itera == 0) {
    for (int k= 0; k < nbElemOptimized; k++)
      ix(k)= D.UI[InitVal_________].D();
  }

  // Use hard coded analytical optimum based on user request for testing purpose
  if (D.UI[ForceSolution___].I() > 0) {
    ix(0)= 6.016;
    ix(1)= 5.309;
    ix(2)= 4.494;
    ix(3)= 3.502;
    ix(4)= 2.153;
  }

  // Set bounds based on user defined move limits and min max values
  for (int k= 0; k < nbElemOptimized; k++) {
    ixmin(k)= std::max(ix(k) - D.UI[MoveLimit_______].D(), D.UI[MinVal__________].D());
    ixmax(k)= std::min(ix(k) + D.UI[MoveLimit_______].D(), D.UI[MaxVal__________].D());
  }

  // Compute objective function and derivative
  double f0val= ix(0) + ix(1) + ix(2) + ix(3) + ix(4);
  for (int k= 0; k < nbElemOptimized; k++)
    df0dx(k)= 1.0;

  // Get value at first iteration for normalization
  static double f0valInit= 1.0;
  if (itera == 0)
    f0valInit= f0val;

  // Normalize objective and derivative
  f0val= f0val / f0valInit - 1.0;
  for (int k= 0; k < nbElemOptimized; k++)
    df0dx(k)= df0dx(k) / f0valInit;

  // Compute constraint function and derivative
  fval(0)= 0.0;
  fval(0)+= 61.0 / std::pow(ix(0), 3.0);
  fval(0)+= 37.0 / std::pow(ix(1), 3.0);
  fval(0)+= 19.0 / std::pow(ix(2), 3.0);
  fval(0)+= 7.0 / std::pow(ix(3), 3.0);
  fval(0)+= 1.0 / std::pow(ix(4), 3.0);
  dfdx(0, 0)= -3.0 * (61.0 / std::pow(ix(0), 4.0));
  dfdx(0, 1)= -3.0 * (37.0 / std::pow(ix(1), 4.0));
  dfdx(0, 2)= -3.0 * (19.0 / std::pow(ix(2), 4.0));
  dfdx(0, 3)= -3.0 * (7.0 / std::pow(ix(3), 4.0));
  dfdx(0, 4)= -3.0 * (1.0 / std::pow(ix(4), 4.0));

  // Normalize constraint and derivative
  fval(0)= fval(0) / D.UI[ConstraintVal___].D() - 1.0;
  for (int k= 0; k < nbElemOptimized; k++)
    dfdx(0, k)= dfdx(0, k) / D.UI[ConstraintVal___].D();

  // Display objective and constraint history in plot area
  if (itera == 0)
    D.Plot.clear();
  D.Plot.resize(2);
  D.Plot[0].name= "Obj";
  D.Plot[1].name= "Const";
  D.Plot[0].val.push_back(f0val);
  D.Plot[1].val.push_back(fval(0));

  // Display solution in scatter area
  if (itera == 0)
    D.Scatter.clear();
  D.Scatter.resize(1);
  D.Scatter[0].name= "Solu";
  D.Scatter[0].val.resize(nbElemOptimized);
  for (int k= 0; k < nbElemOptimized; k++)
    D.Scatter[0].val[k]= std::array<double, 2>({(double)k, ix(k)});

  // Print iteration and current solution values
  printf("itera: %2d ", itera);
  printf("obj: %f ", f0val + 1.0);
  printf("const: %f ", fval(0) + 1.0);
  printf("solu: %f %f %f %f %f\n", ix(0), ix(1), ix(2), ix(3), ix(4));

  // Run one optyimization iteration
  MMA::mmasub(nbConstraints, nbElemOptimized, itera, ix, ixmin, ixmax,
              df0dx, fval, dfdx, xold1, xold2, low, upp, ox);
  itera++;
  ix= ox;
}


// Handle mouse action
void NonLinMMABench::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void NonLinMMABench::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
  if (D.UI[VerboseLevel____].I() >= 5) printf("Animate()\n");
}


// Draw the project
void NonLinMMABench::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;
  if (D.UI[VerboseLevel____].I() >= 5) printf("Draw()\n");
}
