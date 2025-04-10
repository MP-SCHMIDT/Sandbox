#include "NBodyGravDynam.hpp"


// Standard lib
#include <cmath>
#include <numeric>
#include <vector>

// Algo headers
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"

// Project headers
#include "NBodyGravDynam_GPUData.hpp"
#include "NBodyGravDynam_GPUKernel.hpp"


// Link to shared sandbox data
extern Data D;


void NBodyGravDynam::StepSimulation() {

  // Get UI parameters
  const bool lock2D= D.UI[DomainLock2D____].I() == 1;
  const bool torusPos= D.UI[DomainTorusPos__].I() == 1;
  const float dt= D.UI[SimuTimeStep____].F();

  // Optionally sort the bodies
  Timer::PushTimer();
  if (D.UI[SimuBodySort____].I() > 0) {
    MortonSortBodies();
  }
  timerSort= Timer::PopTimer();

  // Semi-implicit Euler integration on CPU
  for (int k= 0; k < D.UI[SimuStepBatch___].I(); k++) {
    ComputeForces(Pos);                         // f(x₀, v₀)
    for (unsigned int k0= 0; k0 < N; k0++) {
      Vel[k0]+= dt * For[k0];                   // v₁ = v₀ + Δt f(x₀, v₀) / m
      Pos[k0]+= dt * Vel[k0];                   // x₁ = x₀ + Δt v₁
      if (lock2D) UtilMake2D(Pos[k0], Vel[k0]); // Apply 2D boundary conditions
      if (torusPos) UtilMakeTorusPos(Pos[k0]);  // Apply periodicity boundary conditions
    }
    simTime+= dt;
  }
}


// Force laws: https://www.desmos.com/calculator/p8ddgpgq5z
void NBodyGravDynam::ComputeForces(const std::vector<Vec::Vec3<float>>& iPos) {
  // Get UI parameters
  const bool torusFor= D.UI[DomainTorusFor__].I() == 1;
  const float gravity= D.UI[SimuTotGravity__].F();
  const float drag= D.UI[SimuDrag________].F();
  const float collision= D.UI[SimuCollision___].F();
  const float distDrag= 8.0f * D.UI[BodyRadius______].F();
  const float distDragInv= 1.0f / distDrag;
  const float distCollision= 2.0f * D.UI[BodyRadius______].F();
  const float distCollisionInv= 1.0f / distCollision;
  const float distCollisionSqr= distCollision * distCollision;
  const float treeTolSqr=  D.UI[SimuTreeTol_____].F() * D.UI[SimuTreeTol_____].F();

  #ifdef TESTING_DISPLAY_FORCES_VECTORS
  ContribPos.clear();
  ContribCount.clear();
  ContribCell.clear();
  #endif

  // Compute forces with N^2 algorithm on CPU
  if (D.UI[SimuMode________].I() == 0) {
    // Sweep through the bodies
    Timer::PushTimer();
    #pragma omp parallel for if (D.UI[SimuMultithread_].I() > 0)
    for (unsigned int k0= 0; k0 < N; k0++) {
      const Vec::Vec3<float> p0(iPos[k0]);
      Vec::Vec3<float> force(0.0f, 0.0f, 0.0f);

      // Sweep through the potential neighbors
      for (unsigned int k1= 0; k1 < N; k1++) {
        if (k0 == k1) continue;

        // Get the neighbor position
        Vec::Vec3<float> p1(iPos[k1]);
        if (p0[0] == p1[0] && p0[1] == p1[1] && p0[2] == p1[2]) continue;
        if (torusFor) UtilMakeTorusFor(p0, p1);
        const Vec::Vec3<float> vec(p1 - p0);
        const float vecNorm= vec.norm();
        const float vecNormInv= 1.0f / vecNorm;
        const float vecNormSqrInv= vecNormInv * vecNormInv;
        const Vec::Vec3<float> vecUnit= vec * vecNormInv;
        
        // Add the gravity force
        force+= vecUnit * gravity * std::min(distCollisionSqr * vecNormSqrInv, vecNorm * distCollisionInv);
        
        // Add the drag force
        if (drag > 0.0f && vecNorm < distDrag)
          force-= (Vel[k0] - Vel[k1]) * drag * (distDrag - vecNorm) * distDragInv;
        
        // Add the collision force
        if (collision > 0.0f && vecNorm < distCollision)
          force-= vecUnit * collision * (distCollision - vecNorm) * distCollisionInv;

        #ifdef TESTING_DISPLAY_FORCES_VECTORS
        if (k0 == N-1) {
          ContribPos.push_back(p1);
          ContribCount.push_back(1);
        }
        #endif
      }

      // Save the resulting force
      For[k0]= force;
    }
    timerForces= Timer::PopTimer();
  }

  // Compute forces with N log(N) algorithm on CPU
  else if (D.UI[SimuMode________].I() == 1) {
    // Build the octree
    Timer::PushTimer();
    BuildTree(iPos);
    timerTree= Timer::PopTimer();
    // Sweep through the bodies
    Timer::PushTimer();
    #pragma omp parallel for if (D.UI[SimuMultithread_].I() > 0)
    for (unsigned int k0= 0; k0 < N; k0++) {
      const Vec::Vec3<float> p0(iPos[k0]);
      Vec::Vec3<float> force(0.0f, 0.0f, 0.0f);

      unsigned int idxCell= 0;
      do {
        // Skip the empty cells
        if (Tree[idxCell].Count == 0) {
          idxCell= Tree[idxCell].Next;
        }
        else {
          const OctreeNode node= Tree[idxCell];
          // Get the cell center of mass
          Vec::Vec3<float> p1(node.AvgPos);
          if (torusFor) UtilMakeTorusFor(p0, p1);
          const Vec::Vec3<float> vec(p1 - p0);

          // Check if the cell can approximate the forces from its sub-tree
          if (node.Child == 0 || node.Count == 1 || treeTolSqr * vec.normSquared() > node.Size * node.Size) {
            if (vec[0] != 0.0f || vec[1] != 0.0f || vec[2] != 0.0f) {
              const float vecNorm= vec.norm();
              const float vecNormInv= 1.0f / vecNorm;
              const float vecNormSqrInv= vecNormInv * vecNormInv;
              const Vec::Vec3<float> vecUnit= vec * vecNormInv;
              const float count= float(node.Count);

              // Add the gravity force
              force+= vecUnit * gravity * count * std::min(distCollisionSqr * vecNormSqrInv, vecNorm * distCollisionInv);

              // Add the drag force
              if (drag > 0.0f && vecNorm < distDrag)
                force-= (Vel[k0] - node.AvgVel) * drag * count * (distDrag - vecNorm) * distDragInv;

              // Add the collision force
              if (collision > 0.0f && vecNorm < distCollision)
                force-= vecUnit * collision * count * (distCollision - vecNorm) * distCollisionInv;

              #ifdef TESTING_DISPLAY_FORCES_VECTORS
              if (k0 == N-1) {
                ContribPos.push_back(p1);
                ContribCount.push_back(node.Count);
                ContribCell.push_back(idxCell);
              }
              #endif
            }
            // The sub-tree has been accounted for, go to the next cell
            idxCell= node.Next;
          }
          else {
            // Go deeper in the sub-tree
            idxCell= node.Child;
          }
        }
      } while (idxCell > 0);

      // Save the resulting force
      For[k0]= force;
    }
    timerForces= Timer::PopTimer();
  }
}


void NBodyGravDynam::SetupGPU() {
  // Check if parameters changed
  if (simTime <= 0.0f)                     OCL_NBody.isMemoryReady= false;
  if (OCL_NBody.N != (int)N)               OCL_NBody.isMemoryReady= false;
  if (D.UI[SimuMode________].hasChanged()) OCL_NBody.isMemoryReady= false;
  if (D.UI[DomainLock2D____].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[DomainTorusPos__].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[DomainTorusFor__].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[SimuTimeStep____].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[SimuTotGravity__].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[SimuDrag________].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[SimuCollision___].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[BodyRadius______].hasChanged()) OCL_NBody.isParamsReady= false;

  // Initialize the device if needed
  if (!OCL_NBody.isDeviceReady) {
    OCL_NBody.isDeviceReady= true;
    OCL_NBody.deviceLink= Device((select_device_with_most_flops(get_devices(false))), NBodyGravDynam_GPUKernel::get_opencl_c_code(), false);
  }

  // Initialize the memory shared between device and host if needed
  if (!OCL_NBody.isMemoryReady) {
    OCL_NBody.isMemoryReady= true;
    OCL_NBody.isParamsReady= false;
    // Allocate the memory shared by host and device
    OCL_NBody.N= N;
    OCL_NBody.arrPos= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
    OCL_NBody.arrVel= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
    OCL_NBody.arrFor= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
    // Set the memory values on the host
    for (int k= 0; k < OCL_NBody.N; k++) {
      OCL_NBody.arrPos[k]= {Pos[k][0], Pos[k][1], Pos[k][2], 0.0f};
      OCL_NBody.arrVel[k]= {Vel[k][0], Vel[k][1], Vel[k][2], 0.0f};
    }
    // Copy the data from host memory to device memory
    OCL_NBody.arrPos.write_to_device();
    OCL_NBody.arrVel.write_to_device();
  }

  // Initialize the kernel parameters
  if (!OCL_NBody.isParamsReady) {
    OCL_NBody.isParamsReady= true;
    // Get the UI parameters
    OCL_NBody.lock2D=        D.UI[DomainLock2D____].I();
    OCL_NBody.torusPos=      D.UI[DomainTorusPos__].I();
    OCL_NBody.torusFor=      D.UI[DomainTorusFor__].I();
    OCL_NBody.timestep=      D.UI[SimuTimeStep____].F();
    OCL_NBody.gravity=       D.UI[SimuTotGravity__].F();
    OCL_NBody.drag=          D.UI[SimuDrag________].F();
    OCL_NBody.collision=     D.UI[SimuCollision___].F();
    OCL_NBody.distDrag=      D.UI[BodyRadius______].F() * 8.0f;
    OCL_NBody.distCollision= D.UI[BodyRadius______].F() * 2.0f;
    // Set the OpenCL kernels
    OCL_NBody.kernelFunc_ComputeForces= Kernel(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, "kernel_ComputeForces",
                                               OCL_NBody.torusFor, OCL_NBody.gravity, OCL_NBody.drag, OCL_NBody.collision,
                                               OCL_NBody.distDrag, OCL_NBody.distCollision,
                                               OCL_NBody.arrPos, OCL_NBody.arrVel, OCL_NBody.arrFor);
    OCL_NBody.kernelFunc_Integrate= Kernel(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, "kernel_Integrate",
                                           OCL_NBody.lock2D, OCL_NBody.torusPos, OCL_NBody.timestep,
                                           OCL_NBody.arrPos, OCL_NBody.arrVel, OCL_NBody.arrFor);
  }
}


void NBodyGravDynam::StepSimulationGPU() {
  // Make sure the GPU is ready for simulation
  SetupGPU();

  // Semi-implicit Euler integration on GPU with optional predictor step
  for (int k= 0; k < std::max(D.UI[SimuStepBatch___].I(), 1); k++) {
    OCL_NBody.kernelFunc_ComputeForces.run();
    OCL_NBody.kernelFunc_Integrate.run();
    simTime+= D.UI[SimuTimeStep____].F();
  }

  // Retrieve data from the GPU device
  OCL_NBody.arrPos.read_from_device();
  OCL_NBody.arrVel.read_from_device();
  OCL_NBody.arrFor.read_from_device();
  for (int k= 0; k < OCL_NBody.N; k++) {
    Pos[k].set(OCL_NBody.arrPos[k].s[0], OCL_NBody.arrPos[k].s[1], OCL_NBody.arrPos[k].s[2]);
    Vel[k].set(OCL_NBody.arrVel[k].s[0], OCL_NBody.arrVel[k].s[1], OCL_NBody.arrVel[k].s[2]);
    For[k].set(OCL_NBody.arrFor[k].s[0], OCL_NBody.arrFor[k].s[1], OCL_NBody.arrFor[k].s[2]);
  }
}
