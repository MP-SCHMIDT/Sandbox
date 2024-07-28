#pragma once

// Standard lib


// Algo headers


// Skeleton Folder
// - Empty template project
class SkeletonFolder
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    AllocTrigger____,
    RefreshTrigger__,
    VerboseLevel____,
  };

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  SkeletonFolder();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void Refresh();
  void ParamChange();
  void KeyPress();
  void MousePress();
  void Animate();
  void Draw();
};
