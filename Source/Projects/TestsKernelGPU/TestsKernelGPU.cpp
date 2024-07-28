#include "TestsKernelGPU.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

// Project headers
#include "TestsKernelGPU_GPUData.hpp"
#include "TestsKernelGPU_GPUKernel.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
TestsKernelGPU::TestsKernelGPU() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void TestsKernelGPU::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ArraySize_______", 1024));
    D.UI.push_back(ParamUI("______________00", NAN));
    D.UI.push_back(ParamUI("ReducSumSize____", 1024));
    D.UI.push_back(ParamUI("______________01", NAN));
    D.UI.push_back(ParamUI("NbParticles_____", 64000));
    D.UI.push_back(ParamUI("InitVel_________", 5.0));
    D.UI.push_back(ParamUI("Timestep________", 0.0004));
    D.UI.push_back(ParamUI("Epsilon_________", 0.001));
    D.UI.push_back(ParamUI("GravCoeff_______", 0.001));
    D.UI.push_back(ParamUI("ColorMode_______", 1));
    D.UI.push_back(ParamUI("ScaleColor______", 0.08));
    D.UI.push_back(ParamUI("ScaleShape______", 0.002));
    D.UI.push_back(ParamUI("______________02", NAN));
    D.UI.push_back(ParamUI("TestParamGPU_00_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_01_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_02_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_03_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_04_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_05_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_06_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_07_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_08_", 0.0));
    D.UI.push_back(ParamUI("TestParamGPU_09_", 0.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));

    D.displayModeLabel[1]= "Point";
    D.displayModeLabel[2]= "Sphere";
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool TestsKernelGPU::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool TestsKernelGPU::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void TestsKernelGPU::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;
}


// Refresh the project
void TestsKernelGPU::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;
}


// Handle UI parameter change
void TestsKernelGPU::ParamChange() {
}


// Handle keypress
void TestsKernelGPU::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  if (D.keyLetterUpperCase == 'A') {
    RunVecAddGPU();
  }

  if (D.keyLetterUpperCase == 'R') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    RunReducSumGPU();
    if (D.UI[VerboseLevel____].I() >= 1) printf("IterT %f\n", Timer::PopTimer());
  }

  if (D.keyLetterUpperCase == 'G') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    StepNBodySimGPU();
    if (D.UI[VerboseLevel____].I() >= 1) printf("IterT %f\n", Timer::PopTimer());
  }

  if (D.keyLetterUpperCase == 'C') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    StepNBodySimCPU();
    if (D.UI[VerboseLevel____].I() >= 1) printf("IterT %f\n", Timer::PopTimer());
  }
}


// Handle mouse action
void TestsKernelGPU::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void TestsKernelGPU::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  StepNBodySimGPU();
}


// Draw the project
void TestsKernelGPU::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Display particles
  if (D.displayMode[1] || D.displayMode[2]) {
    glPointSize(1000.0f * D.UI[ScaleShape______].F());
    if (D.displayMode[1]) glBegin(GL_POINTS);
    else glEnable(GL_LIGHTING);
    for (unsigned int k= 0; k < Pos.size(); k++) {
      float r= 0.5, g= 0.5, b= 0.5;
      if (D.UI[ColorMode_______].I() == 1) {
        Colormap::RatioToJetBrightSmooth(Vel[k].norm() * D.UI[ScaleColor______].F(), r, g, b);
      }
      if (D.UI[ColorMode_______].I() == 2) {
        r= 0.5f + Vel[k][0] * D.UI[ScaleColor______].F();
        g= 0.5f + Vel[k][1] * D.UI[ScaleColor______].F();
        b= 0.5f + Vel[k][2] * D.UI[ScaleColor______].F();
      }
      glColor3f(r, g, b);
      if (D.displayMode[1]) {
        glVertex3fv(Pos[k].array());
      }
      else {
        glPushMatrix();
        glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
        glScalef(D.UI[ScaleShape______].F(), D.UI[ScaleShape______].F(), D.UI[ScaleShape______].F());
        glutSolidSphere(1.0, 16, 8);
        glPopMatrix();
      }
    }
    if (D.displayMode[1]) glEnd();
    else glDisable(GL_LIGHTING);
  }
}


void TestsKernelGPU::RunVecAddGPU() {
  // Initialize the device
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.deviceLink= Device((select_device_with_most_flops(get_devices(D.UI[VerboseLevel____].B()))),
                                TestsKernelGPU_GPUKernel::get_opencl_c_code(), D.UI[VerboseLevel____].B());
  if (D.UI[VerboseLevel____].I() >= 1) printf("CreateDeviceT %f\n", Timer::PopTimer());

  // Get the UI parameters
  OCL_AddVec.N= std::max(D.UI[ArraySize_______].I(), 1);

  // Allocate the memory shared by host and device
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.A= Memory<int>(OCL_AddVec.deviceLink, (ulong)OCL_AddVec.N, 1u, true, true, 0);
  OCL_AddVec.B= Memory<int>(OCL_AddVec.deviceLink, (ulong)OCL_AddVec.N, 1u, true, true, 0);
  OCL_AddVec.C= Memory<int>(OCL_AddVec.deviceLink, (ulong)OCL_AddVec.N, 1u, true, true, 0);
  if (D.UI[VerboseLevel____].I() >= 1) printf("AllocVectorsT %f\n", Timer::PopTimer());

  // Set the OpenCL kernel
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.kernelFuncABC= Kernel(OCL_AddVec.deviceLink, (ulong)OCL_AddVec.N, "kernel_AddVec", OCL_AddVec.A, OCL_AddVec.B, OCL_AddVec.C);
  OCL_AddVec.kernelFuncBCA= Kernel(OCL_AddVec.deviceLink, (ulong)OCL_AddVec.N, "kernel_AddVec", OCL_AddVec.B, OCL_AddVec.C, OCL_AddVec.A);
  if (D.UI[VerboseLevel____].I() >= 1) printf("SetKernelT %f\n", Timer::PopTimer());

  // Set the memory values on the host
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  for (int k= 0; k < OCL_AddVec.N; k++) {
    OCL_AddVec.A[k]= 3;
    OCL_AddVec.B[k]= 2;
  }
  if (D.UI[VerboseLevel____].I() >= 1) printf("FillVectorsT %f\n", Timer::PopTimer());

  // Copy the data from host memory to device memory
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.A.write_to_device();
  OCL_AddVec.B.write_to_device();
  OCL_AddVec.C.write_to_device();
  if (D.UI[VerboseLevel____].I() >= 1) printf("SendVectorsT %f\n", Timer::PopTimer());

  // Run the OpenCL kernel on the data
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.kernelFuncABC.run();
  if (D.UI[VerboseLevel____].I() >= 1) printf("RunKernelT %f\n", Timer::PopTimer());
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.kernelFuncBCA.run();
  if (D.UI[VerboseLevel____].I() >= 1) printf("RunKernelT %f\n", Timer::PopTimer());

  // Copy the data from device memory to host memory
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_AddVec.A.read_from_device();
  OCL_AddVec.B.read_from_device();
  OCL_AddVec.C.read_from_device();
  if (D.UI[VerboseLevel____].I() >= 1) printf("ReceiveVectorsT %f\n", Timer::PopTimer());

  // Print the result
  printf("A[0]= %d B[0]= %d C[0]= %d\n", OCL_AddVec.A[0], OCL_AddVec.B[0], OCL_AddVec.C[0]);
}


void TestsKernelGPU::RunReducSumGPU() {
  OCL_ReducSum.deviceLink= Device((select_device_with_most_flops(get_devices(false))),
                                  TestsKernelGPU_GPUKernel::get_opencl_c_code(), false);
  OCL_ReducSum.N= std::max(D.UI[ReducSumSize____].I(), 1);
  OCL_ReducSum.arr= Memory<int>(OCL_ReducSum.deviceLink, (ulong)OCL_ReducSum.N, 1u, true, true, 0);
  OCL_ReducSum.par= Memory<int>(OCL_ReducSum.deviceLink, OCL_ReducSum.N / WORKGROUP_SIZE + 1u, 1u, true, true, 0);
  OCL_ReducSum.kernelFunc= Kernel(OCL_ReducSum.deviceLink, (ulong)OCL_ReducSum.N, "kernel_ReducSum", OCL_ReducSum.arr, OCL_ReducSum.par);

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  for (int k= 0; k < OCL_ReducSum.N; k++)
    OCL_ReducSum.arr[k]= 1;
  OCL_ReducSum.arr.write_to_device();
  if (D.UI[VerboseLevel____].I() >= 1) printf("ValSetup %f\n", Timer::PopTimer());

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  OCL_ReducSum.kernelFunc.run();
  OCL_ReducSum.par.read_from_device();
  int sumGPU= 0;
  for (int k= 0; k < (int)OCL_ReducSum.par.length(); k++)
    sumGPU+= OCL_ReducSum.par[k];
  if (D.UI[VerboseLevel____].I() >= 1) printf("SumGPUT %f\n", Timer::PopTimer());

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
  int sumCPU= 0;
  for (int k= 0; k < (int)OCL_ReducSum.arr.length(); k++)
    sumCPU+= OCL_ReducSum.arr[k];
  if (D.UI[VerboseLevel____].I() >= 1) printf("SumCPUT %f\n", Timer::PopTimer());

  printf("GPU:%d    CPU:%d\n", sumGPU, sumCPU);
}


void TestsKernelGPU::StepNBodySimGPU() {
  // Check if parameters changed
  if (D.UI[NbParticles_____].hasChanged()) OCL_NBody.isMemoryReady= false;
  if (D.UI[InitVel_________].hasChanged()) OCL_NBody.isMemoryReady= false;
  if (D.UI[Timestep________].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[Epsilon_________].hasChanged()) OCL_NBody.isParamsReady= false;
  if (D.UI[GravCoeff_______].hasChanged()) OCL_NBody.isParamsReady= false;

  // Initialize the device if needed
  if (!OCL_NBody.isDeviceReady) {
    OCL_NBody.isDeviceReady= true;
    OCL_NBody.deviceLink= Device((select_device_with_most_flops(get_devices(false))),
                                 TestsKernelGPU_GPUKernel::get_opencl_c_code(), false);
  }

  // Initialize the memory shared between device and host if needed
  if (!OCL_NBody.isMemoryReady) {
    OCL_NBody.isMemoryReady= true;
    OCL_NBody.isParamsReady= false;
    // Get the UI parameters
    OCL_NBody.N= std::max(D.UI[NbParticles_____].I(), 1);
    // Allocate and initialize bodies with random positions and velocities
    Pos= std::vector<Vec::Vec3<float>>(OCL_NBody.N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    Vel= std::vector<Vec::Vec3<float>>(OCL_NBody.N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    for (int k= 0; k < OCL_NBody.N; k++) {
      Pos[k][0]= (float)D.boxMin[0] + Random::Val(0.0f, 1.0f) * (float)(D.boxMax[0] - D.boxMin[0]);
      Pos[k][1]= (float)D.boxMin[1] + Random::Val(0.0f, 1.0f) * (float)(D.boxMax[1] - D.boxMin[1]);
      Pos[k][2]= (float)D.boxMin[2] + Random::Val(0.0f, 1.0f) * (float)(D.boxMax[2] - D.boxMin[2]);
      Vel[k][0]= std::max(D.UI[InitVel_________].F(), 0.0f) * Random::Val(-1.0f, 1.0f);
      Vel[k][1]= std::max(D.UI[InitVel_________].F(), 0.0f) * Random::Val(-1.0f, 1.0f);
      Vel[k][2]= std::max(D.UI[InitVel_________].F(), 0.0f) * Random::Val(-1.0f, 1.0f);
    }
    // Allocate the memory shared by host and device
    OCL_NBody.PosOld= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
    OCL_NBody.PosNew= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
    OCL_NBody.VelOld= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
    OCL_NBody.VelNew= Memory<cl_float4>(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, 1u, true, true, cl_float4());
  }

  // Initialize the kernel parameters
  if (!OCL_NBody.isParamsReady) {
    OCL_NBody.isParamsReady= true;
    // Get the UI parameters
    OCL_NBody.timestep= std::max(D.UI[Timestep________].F(), 0.0f);
    OCL_NBody.epsilon= std::max(D.UI[Epsilon_________].F(), 0.0f);
    OCL_NBody.gravity= std::max(D.UI[GravCoeff_______].F(), 0.0f);
    // Set the OpenCL kernel
    OCL_NBody.kernelFunc= Kernel(OCL_NBody.deviceLink, (ulong)OCL_NBody.N, "kernel_NBodySim",
                                 OCL_NBody.timestep, OCL_NBody.epsilon, OCL_NBody.gravity,
                                 OCL_NBody.PosOld, OCL_NBody.PosNew, OCL_NBody.VelOld, OCL_NBody.VelNew);
  }

  // Set the memory values on the host
  for (int k= 0; k < OCL_NBody.N; k++) {
    OCL_NBody.PosOld[k]= {Pos[k][0], Pos[k][1], Pos[k][2], 0.0};
    OCL_NBody.VelOld[k]= {Vel[k][0], Vel[k][1], Vel[k][2], 0.0};
  }
  // Copy the data from host memory to device memory
  OCL_NBody.PosOld.write_to_device();
  OCL_NBody.VelOld.write_to_device();
  // Run the OpenCL kernel on the data
  OCL_NBody.kernelFunc.run();
  // Copy the data from device memory to host memory
  OCL_NBody.PosNew.read_from_device();
  OCL_NBody.VelNew.read_from_device();
  // Get the memory values from the host
  for (int k= 0; k < OCL_NBody.N; k++) {
    Pos[k].set(OCL_NBody.PosNew[k].s[0], OCL_NBody.PosNew[k].s[1], OCL_NBody.PosNew[k].s[2]);
    Vel[k].set(OCL_NBody.VelNew[k].s[0], OCL_NBody.VelNew[k].s[1], OCL_NBody.VelNew[k].s[2]);
  }
}


void TestsKernelGPU::StepNBodySimCPU() {
  // Run one iteration of NBody simulation on multithreaded CPU for comparison
  const int N= (int)Pos.size();
  const float dt= D.UI[Timestep________].F();
  const float eps= D.UI[Epsilon_________].F();
  const float grav= D.UI[GravCoeff_______].F();
  std::vector<Vec::Vec3<float>> PosOld= Pos;
  std::vector<Vec::Vec3<float>> VelOld= Vel;
#pragma omp parallel for
  for (int k0= 0; k0 < N; k0++) {          // Loop through bodies
    Vec::Vec3<float> p0= PosOld[k0];       // Get position of current body
    Vec::Vec3<float> v= VelOld[k0];        // Get velocity of current body
    Vec::Vec3<float> a(0.0f, 0.0f, 0.0f);  // Initialize acceleration of current body

    for (int k1= 0; k1 < N; k1++) {                      // Loop through other bodies
      Vec::Vec3<float> p1= PosOld[k1];                   // Get position of other body
      Vec::Vec3<float> rVec= p1 - p0;                    // Get vector from current to other body
      float rNormInv= 1.0f / (rVec.norm() + eps * eps);  // Compute inverse of distance
      a+= grav * rNormInv * rNormInv * rNormInv * rVec;  // Add acceleration toward other body
    }
    v+= dt * a;                        // Integrate velocity
    p0+= dt * v + 0.5f * dt * dt * a;  // Integrate position
    Vel[k0]= v;                        // Save new velocity
    Pos[k0]= p0;                       // Save new position
  }
}
