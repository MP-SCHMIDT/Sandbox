#include "TestsKernelGPU.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// OpenCL lib
#include "OpenCL_Wrapper/opencl.hpp"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Util/Random.hpp"
#include "Util/Timer.hpp"

// Project headers
#include "TestsKernelGPU_Data.hpp"

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
    D.UI.push_back(ParamUI("NbParticles_____", 64000));
    D.UI.push_back(ParamUI("InitVel_________", 5.0));
    D.UI.push_back(ParamUI("Timestep________", 0.0004));
    D.UI.push_back(ParamUI("Epsilon_________", 0.001));
    D.UI.push_back(ParamUI("GravCoeff_______", 0.001));
    D.UI.push_back(ParamUI("ColorMode_______", 1));
    D.UI.push_back(ParamUI("ScaleColor______", 0.08));
    D.UI.push_back(ParamUI("ScaleShape______", 0.002));
    D.UI.push_back(ParamUI("______________01", NAN));
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
  }

  if (D.UI.size() != VerboseLevel____ + 1) {
    printf("[ERROR] Invalid parameter count in UI\n");
  }

  D.boxMin= {0.0, 0.0, 0.0};
  D.boxMax= {1.0, 1.0, 1.0};

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool TestsKernelGPU::CheckAlloc() {
  if (D.UI[NbParticles_____].hasChanged()) isAllocated= false;
  if (D.UI[InitVel_________].hasChanged()) isAllocated= false;
  if (D.UI[Timestep________].hasChanged()) isAllocated= false;
  if (D.UI[Epsilon_________].hasChanged()) isAllocated= false;
  if (D.UI[GravCoeff_______].hasChanged()) isAllocated= false;
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

  // Allocate
  N= std::max(D.UI[NbParticles_____].I(), 1);
  Pos= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  Vel= std::vector<Vec::Vec3<float>>(N, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  isOpenCLReady= false;

  // Set random positions, velocities and types
  for (unsigned int k= 0; k < N; k++) {
    Pos[k][0]= (float)D.boxMin[0] + Random::Val(0.0f, 1.0f) * (float)(D.boxMax[0] - D.boxMin[0]);
    Pos[k][1]= (float)D.boxMin[1] + Random::Val(0.0f, 1.0f) * (float)(D.boxMax[1] - D.boxMin[1]);
    Pos[k][2]= (float)D.boxMin[2] + Random::Val(0.0f, 1.0f) * (float)(D.boxMax[2] - D.boxMin[2]);
    Vel[k][0]= D.UI[InitVel_________].F() * Random::Val(-1.0f, 1.0f);
    Vel[k][1]= D.UI[InitVel_________].F() * Random::Val(-1.0f, 1.0f);
    Vel[k][2]= D.UI[InitVel_________].F() * Random::Val(-1.0f, 1.0f);
  }
}


// Refresh the project
void TestsKernelGPU::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;
}


// Handle keypress
void TestsKernelGPU::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable

  if (key == 'A') {
    // Create the computation device
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    Device_Info device_info(select_device_with_most_flops(get_devices(D.UI[VerboseLevel____].B())));
    Device device(device_info, get_opencl_c_code(), D.UI[VerboseLevel____].B());
    if (D.UI[VerboseLevel____].I() >= 1) printf("CreateDeviceT %f\n", Timer::PopTimer());

    // Allocate the memory shared by host and device
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    const uint N= std::max(D.UI[ArraySize_______].I(), 1);
    Memory<float> A(device, N);
    Memory<float> B(device, N);
    Memory<float> C(device, N);
    if (D.UI[VerboseLevel____].I() >= 1) printf("AllocVectorsT %f\n", Timer::PopTimer());

    // Add the OpenCL kernel
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    Kernel kernel_Add(device, N, "kernel_Add", A, B, C);
    if (D.UI[VerboseLevel____].I() >= 1) printf("SetKernelT %f\n", Timer::PopTimer());

    // Initialize the memory values on the host
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    for (uint n= 0u; n < N; n++) {
      A[n]= 3.0f;
      B[n]= 2.0f;
      C[n]= 1.0f;
    }
    if (D.UI[VerboseLevel____].I() >= 1) printf("FillVectorsT %f\n", Timer::PopTimer());

    // Copy the data from host memory to device memory
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    A.write_to_device();
    B.write_to_device();
    if (D.UI[VerboseLevel____].I() >= 1) printf("SendVectorsT %f\n", Timer::PopTimer());

    // Run the OpenCL kernel on the data
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    kernel_Add.run();
    if (D.UI[VerboseLevel____].I() >= 1) printf("RunKernelT %f\n", Timer::PopTimer());

    // Copy the data from device memory to host memory
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    C.read_from_device();
    if (D.UI[VerboseLevel____].I() >= 1) printf("ReceiveVectorsT %f\n", Timer::PopTimer());

    // Print the resulting value
    if (D.UI[VerboseLevel____].I() >= 1) printf("C[N-1] = %f\n", C[N - 1]);
  }

  if (key == 'C') {
    // Initialize the memory values on the host
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    std::vector<Vec::Vec3<float>> PosOld= Pos;
    std::vector<Vec::Vec3<float>> VelOld= Vel;
    if (D.UI[VerboseLevel____].I() >= 1) printf("FillVectorsT %f\n", Timer::PopTimer());

    // Run the OpenCL kernel on the data
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    float dt= D.UI[Timestep________].F();
    float eps= D.UI[Epsilon_________].F();
    float grav= D.UI[GravCoeff_______].F();
    for (unsigned int k0= 0; k0 < N; k0++) {  // Loop through bodies
      Vec::Vec3<float> p0= PosOld[k0];        // Get position of current body
      Vec::Vec3<float> v= VelOld[k0];         // Get velocity of current body
      Vec::Vec3<float> a(0.0f, 0.0f, 0.0f);   // Initialize acceleration of current body

      for (unsigned int k1= 0; k1 < N; k1++) {             // Loop through other bodies
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
    if (D.UI[VerboseLevel____].I() >= 1) printf("RunKernelT %f\n", Timer::PopTimer());
  }

  // https://dournac.org/info/nbody_tutorial
  if (key == 'G') {
    if (!isOpenCLReady) {
      isOpenCLReady= true;
      // Create the computation device
      if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
      device= Device((select_device_with_most_flops(get_devices(D.UI[VerboseLevel____].B()))), get_opencl_c_code(), D.UI[VerboseLevel____].B());
      if (D.UI[VerboseLevel____].I() >= 1) printf("CreateDeviceT %f\n", Timer::PopTimer());

      // Allocate the memory shared by host and device
      if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
      PosOld= Memory<cl_float4>(device, N, 1u, true, true, cl_float4());
      PosNew= Memory<cl_float4>(device, N, 1u, true, true, cl_float4());
      VelOld= Memory<cl_float4>(device, N, 1u, true, true, cl_float4());
      VelNew= Memory<cl_float4>(device, N, 1u, true, true, cl_float4());
      if (D.UI[VerboseLevel____].I() >= 1) printf("AllocVectorsT %f\n", Timer::PopTimer());

      // Add the OpenCL kernel
      if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
      kernel_NBody= Kernel(device, N, "kernel_NBody",
                           N, D.UI[Timestep________].F(), D.UI[Epsilon_________].F(), D.UI[GravCoeff_______].F(),
                           PosOld, PosNew, VelOld, VelNew);
      if (D.UI[VerboseLevel____].I() >= 1) printf("SetKernelT %f\n", Timer::PopTimer());
    }

    // Initialize the memory values on the host
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    for (uint n= 0u; n < N; n++) {
      PosOld[n]= {Pos[n][0], Pos[n][1], Pos[n][2], 0.0};
      VelOld[n]= {Vel[n][0], Vel[n][1], Vel[n][2], 0.0};
    }
    if (D.UI[VerboseLevel____].I() >= 1) printf("FillVectorsT %f\n", Timer::PopTimer());

    // Copy the data from host memory to device memory
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    PosOld.write_to_device();
    VelOld.write_to_device();
    if (D.UI[VerboseLevel____].I() >= 1) printf("SendVectorsT %f\n", Timer::PopTimer());

    // Run the OpenCL kernel on the data
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    kernel_NBody.run();
    if (D.UI[VerboseLevel____].I() >= 1) printf("RunKernelT %f\n", Timer::PopTimer());

    // Copy the data from device memory to host memory
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    PosNew.read_from_device();
    VelNew.read_from_device();
    if (D.UI[VerboseLevel____].I() >= 1) printf("ReceiveVectorsT %f\n", Timer::PopTimer());

    // Read the memory values on the host
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    for (uint n= 0u; n < N; n++) {
      Pos[n].set(PosNew[n].s[0], PosNew[n].s[1], PosNew[n].s[2]);
      Vel[n].set(VelNew[n].s[0], VelNew[n].s[1], VelNew[n].s[2]);
    }
    if (D.UI[VerboseLevel____].I() >= 1) printf("ReadVectorsT %f\n", Timer::PopTimer());
  }
}


// Handle mouse action
void TestsKernelGPU::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
}


// Animate the project
void TestsKernelGPU::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
}


// Draw the project
void TestsKernelGPU::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Display particles

  if (D.displayMode1) {
    glPointSize(1000.0f * D.UI[ScaleShape______].F());
    glBegin(GL_POINTS);
    for (unsigned int k= 0; k < N; k++) {
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
      glVertex3fv(Pos[k].array());
    }
    glEnd();
  }

  if (!D.displayMode2) {
    glEnable(GL_LIGHTING);
    for (unsigned int k= 0; k < N; k++) {
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
      glPushMatrix();
      glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
      glScalef(D.UI[ScaleShape______].F(), D.UI[ScaleShape______].F(), D.UI[ScaleShape______].F());
      glutSolidSphere(1.0, 16, 8);
      glPopMatrix();
    }
    glDisable(GL_LIGHTING);
  }
}
