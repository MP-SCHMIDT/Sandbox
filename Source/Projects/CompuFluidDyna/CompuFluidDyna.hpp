#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Fluid simulation code
// - Eulerian voxel grid
// - Handles 1D, 2D and 3D transparently
// - Linear solve in implicit diffusion step for viscosity and smoke spread/mixing
// - Linear solve in implicit pressure computation and projection to enforce mass conservation
// - Solves all linear systems with a custom matrixless diagonal preconditioned conjugate gradient
// - Semi Lagrangian backtracing for velocity and smoke advection
// - Uses iterative MackCormack backtracking scheme to achieve 2nd order accuracy in advection steps
// - Reinjects dissipated vorticity at smallest scale using vorticity confinement approach
// - Handles arbitrary boundary conditions and obstacles in the simulation domain using boolean flag fields
// - Validated on Re < 1000 in lid-driven cavity flow, Poiseuille, Couette and venturi benchmarks
// - Uses SI units
//
// References for fluid flow photographs, scenarios and visual comparison
// http://courses.washington.edu/me431/handouts/Album-Fluid-Motion-Van-Dyke.pdf
//
// References for "stable fluid" method
// http://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
// https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf
// https://www.dgp.toronto.edu/public_user/stam/reality/Research/pub.html
// https://fr.wikipedia.org/wiki/Stable-Fluids
// http://www.dgp.utoronto.ca/~stam/reality/Talks/FluidsTalk/FluidsTalkNotes.pdf
// https://www.youtube.com/watch?v=qsYE1wMEMPA theory simple explanation
// https://www.youtube.com/watch?v=iKAVRgIrUOU JS, Matthias Müller, slightly different approach for pressure
// https://www.youtube.com/watch?v=wbYe58NGJJI python
// https://github.com/NiallHornFX/StableFluids3D-GL/blob/master/src/fluidsolver3d.cpp
//
// TODO Test heuristic optimization of solid regions
// Optimisation of pipes with constant diameter using the heuristic optimality criterion - David Blacher, Michael Harasek
// https://open-research-europe.ec.europa.eu/articles/3-156
// https://github.com/AIT-LKR/SEROS
//
// TODO switch to continusous Solid field and implement Brinkman penalization ?
// On the Calculation of the Brinkman Penalization Term in Density-Based Topology Optimization of Fluid-Dependent Problems - Mohamed Abdelhamid, Aleksander Czekanski
// https://arxiv.org/pdf/2302.14156.pdf
// Towards improved porous models for solid/fluid topology optimization - Maarten J. B. Theulings, Matthijs Langelaar, Fred van Keulen & Robert Maas
// https://link.springer.com/article/10.1007/s00158-023-03570-4
//
// TODO Evaluate feasibility of adjoint analysis to get gradients for optim
// Topology optimization of unsteady flow problems using the lattice Boltzmann method - Sebastian Arlund Nørgaard, Ole Sigmund, Boyan S Lazarov
// https://www.researchgate.net/publication/287621818_Topology_optimization_of_unsteady_flow_problems_using_the_lattice_Boltzmann_method
// https://github.com/northgaard/LBM-TopOpt
// A detailed introduction to density‑based topology optimisation of fluid flow problems with implementation in MATLAB - Joe Alexandersen
// https://link.springer.com/article/10.1007/s00158-022-03420-9
// https://github.com/sdu-multiphysics/topflow
class CompuFluidDyna
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ResolutionX_____,
    ResolutionY_____,
    ResolutionZ_____,
    VoxelSize_______,
    Scenario________,
    InputFile_______,
    TimeStep________,
    Multithread_____,
    SolvMaxIter_____,
    SolvType________,
    SolvTolRhs______,
    SolvTolRel______,
    SolvTolAbs______,
    CoeffGravi______,
    CoeffAdvec______,
    CoeffAdvecTol___,
    CoeffDiffuS_____,
    CoeffDiffuV_____,
    CoeffVorti______,
    CoeffProj_______,
    PorosSmoothIt___,
    PorosMinThresh__,
    PorosOffset_____,
    DarcyMinResist__,
    DarcyMaxResist__,
    DarcyPenal______,
    BCVelX__________,
    BCVelY__________,
    BCVelZ__________,
    BCPres__________,
    BCSmok__________,
    BCSmokTime______,
    ObjectPosX______,
    ObjectPosY______,
    ObjectPosZ______,
    ObjectSize0_____,
    ObjectSize1_____,
    ParticCount_____,
    ParticDuration__,
    ScaleFactor_____,
    ColorFactor_____,
    ColorThresh_____,
    ColorMode_______,
    SliceDim________,
    SlicePlotX______,
    SlicePlotY______,
    SlicePlotZ______,
    PlotSolve_______,
    PlotMFR0Offset__,
    PlotMFR1Offset__,
    VerboseLevel____,
  };

  // List of field IDs
  enum FieldID
  {
    IDSmok,
    IDVelX,
    IDVelY,
    IDVelZ,
    IDPres,
  };

  // Problem dimensions
  int nX;
  int nY;
  int nZ;
  float voxSize;
  float simTime;

  // Fluid properties
  float fluidDensity;

  // Fields for scenario setup
  std::vector<std::vector<std::vector<bool>>> Solid;
  std::vector<std::vector<std::vector<bool>>> VelBC;
  std::vector<std::vector<std::vector<bool>>> PreBC;
  std::vector<std::vector<std::vector<bool>>> SmoBC;
  std::vector<std::vector<std::vector<float>>> VelXForced;
  std::vector<std::vector<std::vector<float>>> VelYForced;
  std::vector<std::vector<std::vector<float>>> VelZForced;
  std::vector<std::vector<std::vector<float>>> PresForced;
  std::vector<std::vector<std::vector<float>>> SmokForced;
  std::vector<std::vector<std::vector<float>>> Poros;

  // Fields for scenario run
  std::vector<std::vector<std::vector<float>>> Dum0;
  std::vector<std::vector<std::vector<float>>> Dum1;
  std::vector<std::vector<std::vector<float>>> Dum2;
  std::vector<std::vector<std::vector<float>>> Dum3;
  std::vector<std::vector<std::vector<float>>> Dum4;
  std::vector<std::vector<std::vector<float>>> Vort;
  std::vector<std::vector<std::vector<float>>> Pres;
  std::vector<std::vector<std::vector<float>>> Dive;
  std::vector<std::vector<std::vector<float>>> Smok;
  std::vector<std::vector<std::vector<float>>> VelX;
  std::vector<std::vector<std::vector<float>>> VelY;
  std::vector<std::vector<std::vector<float>>> VelZ;
  std::vector<std::vector<std::vector<float>>> CurX;
  std::vector<std::vector<std::vector<float>>> CurY;
  std::vector<std::vector<std::vector<float>>> CurZ;
  std::vector<std::vector<std::vector<float>>> AdvX;
  std::vector<std::vector<std::vector<float>>> AdvY;
  std::vector<std::vector<std::vector<float>>> AdvZ;

  // Array for streamline or particle display
  std::vector<Vec::Vec3<float>> ParticlesPos;
  std::vector<float> ParticlesAge;

  // Initialization and display functions
  void InitializeScenario();
  void UpdateUIData();

  // Simulation functions
  void RunSimulationStep();
  void UpdateSolidFieldFromPoros();
  void ApplyBC(const int iFieldID, std::vector<std::vector<std::vector<float>>>& ioField);
  void ExternalGravityForce();
  void ExternalDarcyForce();
  void ProjectField(const int iMaxIter, const float iTimeStep,
                    std::vector<std::vector<std::vector<float>>>& ioVelX,
                    std::vector<std::vector<std::vector<float>>>& ioVelY,
                    std::vector<std::vector<std::vector<float>>>& ioVelZ);
  float TrilinearInterpolation(const float iPosX, const float iPosY, const float iPosZ,
                               const std::vector<std::vector<std::vector<float>>>& iFieldRef);
  void AdvectField(const int iFieldID, const float iTimeStep,
                   const std::vector<std::vector<std::vector<float>>>& iVelX,
                   const std::vector<std::vector<std::vector<float>>>& iVelY,
                   const std::vector<std::vector<std::vector<float>>>& iVelZ,
                   std::vector<std::vector<std::vector<float>>>& ioField);
  void VorticityConfinement(const float iTimeStep, const float iVortiCoeff,
                            std::vector<std::vector<std::vector<float>>>& ioVelX,
                            std::vector<std::vector<std::vector<float>>>& ioVelY,
                            std::vector<std::vector<std::vector<float>>>& ioVelZ);
  void ComputeParticlesMovement();

  // Vector field operators
  void ComputeVectorFieldDivergence(const std::vector<std::vector<std::vector<float>>>& iVecX,
                                    const std::vector<std::vector<std::vector<float>>>& iVecY,
                                    const std::vector<std::vector<std::vector<float>>>& iVecZ,
                                    std::vector<std::vector<std::vector<float>>>& oDiv);
  void ComputeVectorFieldCurl(const std::vector<std::vector<std::vector<float>>>& iVecX,
                              const std::vector<std::vector<std::vector<float>>>& iVecY,
                              const std::vector<std::vector<std::vector<float>>>& iVecZ,
                              std::vector<std::vector<std::vector<float>>>& oCurlX,
                              std::vector<std::vector<std::vector<float>>>& oCurlY,
                              std::vector<std::vector<std::vector<float>>>& oCurlZ);
  void ComputeVectorFieldNorm(const std::vector<std::vector<std::vector<float>>>& iVecX,
                              const std::vector<std::vector<std::vector<float>>>& iVecY,
                              const std::vector<std::vector<std::vector<float>>>& iVecZ,
                              std::vector<std::vector<std::vector<float>>>& oNorm);

  // Linear solver function helpers
  void ImplicitFieldAdd(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                        const std::vector<std::vector<std::vector<float>>>& iFieldB,
                        std::vector<std::vector<std::vector<float>>>& oField);
  void ImplicitFieldSub(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                        const std::vector<std::vector<std::vector<float>>>& iFieldB,
                        std::vector<std::vector<std::vector<float>>>& oField);
  void ImplicitFieldScale(const float iVal,
                          const std::vector<std::vector<std::vector<float>>>& iField,
                          std::vector<std::vector<std::vector<float>>>& oField);
  float ImplicitFieldDotProd(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                             const std::vector<std::vector<std::vector<float>>>& iFieldB);
  void ImplicitFieldLaplacianMatMult(const int iFieldID, const float iTimeStep, const float iDiffuCoeff, const bool iPrecondMode,
                                     const std::vector<std::vector<std::vector<float>>>& iField,
                                     std::vector<std::vector<std::vector<float>>>& oField);

  // Linear solver functions
  void LinearSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                   const std::vector<std::vector<std::vector<float>>>& iField,
                   std::vector<std::vector<std::vector<float>>>& ioField);
  void PreconditionedConjugateGradientSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                                            const std::vector<std::vector<std::vector<float>>>& iField,
                                            std::vector<std::vector<std::vector<float>>>& ioField);
  void ConjugateGradientSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                              const std::vector<std::vector<std::vector<float>>>& iField,
                              std::vector<std::vector<std::vector<float>>>& ioField);
  void GradientDescentSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                            const std::vector<std::vector<std::vector<float>>>& iField,
                            std::vector<std::vector<std::vector<float>>>& ioField);
  void GaussSeidelSolve(const int iFieldID, const int iMaxIter, const float iTimeStep, const float iDiffuCoeff,
                        const std::vector<std::vector<std::vector<float>>>& iField,
                        std::vector<std::vector<std::vector<float>>>& ioField);

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  CompuFluidDyna();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
