#include "SkeletonFolder.hpp"


// Standard lib

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
SkeletonFolder::SkeletonFolder() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void SkeletonFolder::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("AllocTrigger____", 0.0));
    D.UI.push_back(ParamUI("RefreshTrigger__", 0.0));
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
bool SkeletonFolder::CheckAlloc() {
  if (D.UI[AllocTrigger____].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool SkeletonFolder::CheckRefresh() {
  if (D.UI[RefreshTrigger__].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void SkeletonFolder::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Allocate()\n");
}


// Refresh the project
void SkeletonFolder::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Refresh()\n");
}


// Handle keypress
void SkeletonFolder::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable

  if (D.UI[VerboseLevel____].I() >= 5) printf("KeyPress()\n");
}


// Handle mouse action
void SkeletonFolder::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
}


// Animate the project
void SkeletonFolder::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
  if (D.UI[VerboseLevel____].I() >= 5) printf("Animate()\n");
}


// Draw the project
void SkeletonFolder::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;
  if (D.UI[VerboseLevel____].I() >= 5) printf("Draw()\n");
}
