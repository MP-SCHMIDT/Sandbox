#pragma once


class ImageExtruMesh
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    ResolutionX_____,
    ResolutionY_____,
    ResolutionZ_____,
    SizeX___________,
    SizeY___________,
    SizeZ___________,
    CenterX_________,
    CenterY_________,
    CenterZ_________,
    BaseRelHeight___,
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
