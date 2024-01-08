#pragma once


class ImageExtruMesh
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    InputFile_______,
    DomainW_________,
    DomainH_________,
    SmoothIter______,
    Isovalue________,
    VerboseLevel____,
  };

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  ImageExtruMesh();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void KeyPress(const unsigned char key);
  void Refresh();
  void Animate();
  void Draw();
};
