#include "NBodyGravDynam.hpp"


// Standard lib
#include <format>
#include <random>

// Algo headers
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
NBodyGravDynam::NBodyGravDynam() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void NBodyGravDynam::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("BodyCount_______", 20000)); // Number of bodies in the simulation
    D.UI.push_back(ParamUI("BodyRandSeed____", 1));     // Fixed seed for the random number generators
    D.UI.push_back(ParamUI("BodyInitLayout__", 1));     // Initial layout of bodies according to different scenarios
    D.UI.push_back(ParamUI("BodyInitVel_____", 10));    // Reference initial velocity
    D.UI.push_back(ParamUI("BodyRadius______", 0.003)); // Body radius used for collision and drag
    D.UI.push_back(ParamUI("______________00", NAN));   //
    D.UI.push_back(ParamUI("DomainLock2D____", 0));     // Constrain the simulation to 2D
    D.UI.push_back(ParamUI("DomainTorusPos__", 0));     // Make positions follow a periodic unit box domain
    D.UI.push_back(ParamUI("DomainTorusFor__", 0));     // Make forces follow a periodic unit box domain
    D.UI.push_back(ParamUI("______________01", NAN));   //
    D.UI.push_back(ParamUI("TreeMaxDepth____", 16));    // Maximum tree depth
    D.UI.push_back(ParamUI("TreeInfiniteBox_", 1));     // Automaticaly compute and fit the tree to the cloud bounding box
    D.UI.push_back(ParamUI("______________02", NAN));   //
    D.UI.push_back(ParamUI("SimuMode________", 1));     // Simulation mode: 0= N^2 on CPU, 1= N log(N) on CPU, 2= N^2 on GPU
    D.UI.push_back(ParamUI("SimuTreeTol_____", 0.7));   // Tree cell distance tolerance (0.0 becomes N^2 algorithm, 0.5 is the original Barnes-Hut value)
    D.UI.push_back(ParamUI("SimuTimeStep____", 0.001)); // Simulation timestep
    D.UI.push_back(ParamUI("SimuStepBatch___", 1));     // Number of simultation step computed between each draw
    D.UI.push_back(ParamUI("SimuTotGravity__", 1.0));   // Total gravity force exerted by the cloud
    D.UI.push_back(ParamUI("SimuDrag________", 0.1));   // Pairwise drag force exerted between the bodies
    D.UI.push_back(ParamUI("SimuCollision___", 128.0)); // Pairwise collision force exerted between the bodies
    D.UI.push_back(ParamUI("SimuBodySort____", 1));     // Sort the bodies along Morton curve to accelerate tree traversal
    D.UI.push_back(ParamUI("SimuMultithread_", 1));     // Enable multithreading for force calculation
    D.UI.push_back(ParamUI("______________03", NAN));   //
    D.UI.push_back(ParamUI("ColorMode_______", 1));     //
    D.UI.push_back(ParamUI("ColorFactor_____", 0.1));   //
    D.UI.push_back(ParamUI("ScaleFactor_____", 0.1));   //
    D.UI.push_back(ParamUI("SphereSimple____", 2));     //
    D.UI.push_back(ParamUI("ShowEmptyCells__", 0));     //
    D.UI.push_back(ParamUI("______________04", NAN));   //
    D.UI.push_back(ParamUI("TestParamNBS_00_", 0.0));   //
    D.UI.push_back(ParamUI("TestParamNBS_01_", 0.0));   //
    D.UI.push_back(ParamUI("TestParamNBS_02_", 0.0));   //
    D.UI.push_back(ParamUI("TestParamNBS_03_", 0.0));   //
    D.UI.push_back(ParamUI("TestParamNBS_04_", 0.0));   //
    D.UI.push_back(ParamUI("VerboseLevel____", 0));     //

    D.displayModeLabel[1]= "Bodies Pos";
    D.displayModeLabel[2]= "Bodies Vel";
    D.displayModeLabel[3]= "Bodies Order";
    D.displayModeLabel[4]= "Octree";
    D.displayModeLabel[5]= "Octree AvgPos";
    D.displayModeLabel[6]= "Octree AvgVel";
    D.displayModeLabel[7]= "Octree Order";
    #ifdef TESTING_DISPLAY_FORCES_VECTORS
    D.displayModeLabel[8]= "Approx Force";
    D.displayModeLabel[9]= "Approx Source";
    #endif
    D.displayMode[2]= false;
    D.displayMode[3]= false;
    D.displayMode[4]= false;
    D.displayMode[5]= false;
    D.displayMode[6]= false;
    D.displayMode[7]= false;
    D.displayMode[8]= false;
    D.displayMode[9]= false;
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool NBodyGravDynam::CheckAlloc() {
  if (D.UI[BodyCount_______].hasChanged()) isAllocated= false;
  if (D.UI[BodyRandSeed____].hasChanged()) isAllocated= false;
  if (D.UI[BodyInitLayout__].hasChanged()) isAllocated= false;
  if (D.UI[BodyInitVel_____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool NBodyGravDynam::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void NBodyGravDynam::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Set the box domain
  D.boxMin= {0.0f, 0.0f, 0.0f};
  D.boxMax= {1.0f, 1.0f, 1.0f};

  // Get UI parameters
  N= std::max(D.UI[BodyCount_______].I(), 1);
  simTime= 0.0f;

  // Allocate data
  Pos= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  Vel= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  For= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));

  // Initialize random bodies with fixed seed
  if (D.UI[BodyInitLayout__].I() == 0) {
    std::default_random_engine rng(D.UI[BodyRandSeed____].I());
    std::uniform_real_distribution<float> uDist(0.0f, 1.0f);
    std::uniform_real_distribution<float> uDistSym(-1.0f, 1.0f);
    for (unsigned int k0= 0; k0 < N; k0++) {
      for (unsigned int dim= 0; dim < 3; dim++) {
        Pos[k0][dim]= uDist(rng);
        Vel[k0][dim]= uDistSym(rng) * D.UI[BodyInitVel_____].F();
      }
    }
  }
  else if (D.UI[BodyInitLayout__].I() == 1) {
    std::default_random_engine rng(D.UI[BodyRandSeed____].I());
    std::uniform_real_distribution<float> uDist(0.0f, 1.0f);
    std::uniform_real_distribution<float> uDistSym(-1.0f, 1.0f);
    std::normal_distribution<float> nDist;
    const Vec::Vec3<float> center(0.5f, 0.5f, 0.5f);
    const Vec::Vec3<float> spinAxis(1.0f, 0.0f, 0.0f);
    for (unsigned int k0= 0; k0 < N; k0++) {
      for (unsigned int dim= 0; dim < 3; dim++)
        Pos[k0][dim]= nDist(rng);
      Pos[k0]= center + 0.25f * std::cbrt(uDist(rng)) * Pos[k0].normalized();
      const Vec::Vec3<float> tangent= spinAxis.cross(Pos[k0] - center) / spinAxis.cross(Pos[k0] - center).norm();
      Vel[k0]= tangent * D.UI[BodyInitVel_____].F() * (0.5f + 0.2f * uDistSym(rng)) * (center - Pos[k0]).norm();
    }
  }
  else if (D.UI[BodyInitLayout__].I() == 2) {
    std::default_random_engine rng(D.UI[BodyRandSeed____].I());
    std::uniform_real_distribution<float> uDist(0.0f, 1.0f);
    std::uniform_real_distribution<float> uDistSym(-1.0f, 1.0f);
    std::normal_distribution<float> nDist;
    const Vec::Vec3<float> center(0.5f, 0.5f, 0.5f);
    const Vec::Vec3<float> spinAxis(1.0f, 0.0f, 0.0f);
    for (unsigned int k0= 0; k0 < N; k0++) {
      for (unsigned int dim= 0; dim < 3; dim++)
        Pos[k0][dim]= center[dim] + 0.1f * nDist(rng) * nDist(rng);
      const Vec::Vec3<float> tangent= spinAxis.cross(Pos[k0] - center) / spinAxis.cross(Pos[k0] - center).norm();
      Vel[k0]= tangent * D.UI[BodyInitVel_____].F() * (1.0f + 0.1f * nDist(rng)) * std::sqrt(std::max(D.UI[SimuTotGravity__].F() - (center - Pos[k0]).norm(), 0.0f));
    }
  }

  // Apply boundary conditions
  for (unsigned int k0= 0; k0 < N; k0++) {
    if (D.UI[DomainLock2D____].I() == 1) UtilMake2D(Pos[k0], Vel[k0]);
    if (D.UI[DomainTorusPos__].I() == 1) UtilMakeTorusPos(Pos[k0]);
  }
}


// Refresh the project
void NBodyGravDynam::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Run one iteration to set up the data and visualization
  Animate();
}


// Handle UI parameter change
void NBodyGravDynam::ParamChange() {
}


// Handle keypress
void NBodyGravDynam::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Handle mouse action
void NBodyGravDynam::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  Vec::Vec3<float> mousePos((float)D.mouseProjX[0], (float)D.mouseProjX[1], (float)D.mouseProjX[2]);

  if (D.mouseMiddleButtonState == 1) {
    N++;
    Pos.push_back(mousePos);
    Vel.push_back({0.0f, 0.0f, 0.0f});
    For.push_back({0.0f, 0.0f, 0.0f});
  }
  else if (D.mouseMiddleButtonState == 2) {
    Vel[N-1]= (mousePos - Pos[N-1]) / D.UI[ScaleFactor_____].F();
  }

}


// Animate the project
void NBodyGravDynam::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  Timer::PushTimer();
  if      (D.UI[SimuMode________].I() == 0) StepSimulation();
  else if (D.UI[SimuMode________].I() == 1) StepSimulation();
  else                                      StepSimulationGPU();
  timerSimu= Timer::PopTimer();

  if (D.UI[VerboseLevel____].I() >= 1) printf("\n");
}


// Draw the project
void NBodyGravDynam::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  Timer::PushTimer();
  DrawScene();
  timerDraw= Timer::PopTimer();

  D.Status.clear();
  D.Status.push_back(std::format("SimTime:{:.3f}s", simTime));
  D.Status.push_back(std::format("Bodies:{}={}MB", N, N * 3 * sizeof(Vec::Vec3<float>) / 1000000));
  D.Status.push_back(std::format("Tree cells: {}={}MB", Tree.size(), Tree.size() * sizeof(NBodyGravDynam::OctreeNode) / 1000000));
  D.Status.push_back(std::format("TSort: {:.3f}s", timerSort));
  D.Status.push_back(std::format("TTree: {:.3f}s", timerTree));
  D.Status.push_back(std::format("TForces: {:.3f}s", timerForces));
  D.Status.push_back(std::format("TSim: {:.3f}s", timerSimu));
  D.Status.push_back(std::format("TDraw: {:.3f}s", timerDraw));
}
