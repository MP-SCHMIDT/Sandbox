#include "TestsKernelGPU.hpp"


// Standard lib
#include <cmath>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// OpenCL lib
#include "OpenCL_Wrapper/opencl.hpp"

// Algo headers
#include "Util/Timer.hpp"

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


// Handle keypress
void TestsKernelGPU::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable

  if (key == 'T') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

    Device_Info device_info= select_device_with_most_flops(get_devices(D.UI[VerboseLevel____].B()));  // compile OpenCL C code for the fastest available device
    Device device(device_info, get_opencl_c_code(), D.UI[VerboseLevel____].B());
    if (D.UI[VerboseLevel____].I() >= 1) printf("CreateDeviceT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    const uint N= D.UI[ArraySize_______].I();  // size of vectors
    Memory<float> A(device, N);                // allocate memory on both host and device
    Memory<float> B(device, N);
    Memory<float> C(device, N);
    if (D.UI[VerboseLevel____].I() >= 1) printf("AllocVectorsT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    Kernel add_kernel(device, N, "add_kernel", A, B, C);  // kernel that runs on the device
    if (D.UI[VerboseLevel____].I() >= 1) printf("AddKernelT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    for (uint n= 0u; n < N; n++) {
      A[n]= 3.0f;  // initialize memory
      B[n]= 2.0f;
      C[n]= 1.0f;
    }
    if (D.UI[VerboseLevel____].I() >= 1) printf("FillVectorsT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) printf("C[N-1] = %f\n", C[N - 1]);

    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    A.write_to_device();  // copy data from host memory to device memory
    B.write_to_device();
    if (D.UI[VerboseLevel____].I() >= 1) printf("SendVectorsT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    add_kernel.run();  // run add_kernel on the device
    if (D.UI[VerboseLevel____].I() >= 1) printf("RunKernelT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    C.read_from_device();  // copy data from device memory to host memory
    if (D.UI[VerboseLevel____].I() >= 1) printf("ReceiveVectorsT %f\n", Timer::PopTimer());

    if (D.UI[VerboseLevel____].I() >= 1) printf("C[N-1] = %f\n", C[N - 1]);
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
}
