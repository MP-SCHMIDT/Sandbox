#pragma once

// Standard lib


// Algo headers


// Skeleton Folder
// - Empty template project
class NonLinMMABench
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ResetOptim______,
    ForceSolution___,
    InitVal_________,
    MinVal__________,
    MaxVal__________,
    MoveLimit_______,
    ConstraintVal___,
    VerboseLevel____,
  };

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  NonLinMMABench();

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
