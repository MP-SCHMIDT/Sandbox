#include "ImageExtruMesh.hpp"


// Standard lib
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "FileIO/FileInput.hpp"
#include "FileIO/FileOutput.hpp"
#include "Geom/MarchingCubes.hpp"
#include "Math/Field.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
ImageExtruMesh::ImageExtruMesh() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void ImageExtruMesh::SetActiveProject() {
  if (!isActivProj) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ResolutionX_____", 100));
    D.UI.push_back(ParamUI("ResolutionY_____", 100));
    D.UI.push_back(ParamUI("SizeX___________", 10.0));
    D.UI.push_back(ParamUI("SizeY___________", 10.0));
    D.UI.push_back(ParamUI("SizeZ___________", 1.0));
    D.UI.push_back(ParamUI("SmoothIter______", 3));
    D.UI.push_back(ParamUI("Isovalue________", 0.5));
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
bool ImageExtruMesh::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool ImageExtruMesh::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void ImageExtruMesh::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;
  if (D.UI[VerboseLevel____].I() >= 5) printf("ImageExtruMesh::Allocate()\n");
}


// Refresh the project
void ImageExtruMesh::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;
  if (D.UI[VerboseLevel____].I() >= 5) printf("ImageExtruMesh::Refresh()\n");

  // Set the bounding box
  D.boxMin= {-0.5 * D.UI[SizeX___________].D(), -0.5 * D.UI[SizeY___________].D(), -0.5 * D.UI[SizeZ___________].D()};
  D.boxMax= {0.5 * D.UI[SizeX___________].D(), 0.5 * D.UI[SizeY___________].D(), 0.5 * D.UI[SizeZ___________].D()};

  // Get and check dimensions
  const int nX= D.UI[ResolutionX_____].I();
  const int nY= D.UI[ResolutionY_____].I();
  if (nX < 1 || nY < 1) return;

  // Load bitmap file
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  FileInput::LoadImageBMPFile("./FileInput/Images/Logo.bmp", imageRGBA, false);

  // Get bitmap file dimensions
  int nW, nH;
  Field::GetFieldDimensions(imageRGBA, nW, nH);
  if (nW < 1 || nH < 1) return;

  // Sample field from image red channel with zero order interpolation (nearest neighbor) scheme
  std::vector<std::vector<std::vector<double>>> field= Field::AllocField3D(nX, nY, 1, 0.0);
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      const float posW= (float)(nW - 1) * ((float)x + 0.5f) / (float)nX;
      const float posH= (float)(nH - 1) * ((float)y + 0.5f) / (float)nY;
      const int idxPixelW= std::min(std::max((int)std::round(posW), 0), nW - 1);
      const int idxPixelH= std::min(std::max((int)std::round(posH), 0), nH - 1);
      field[x][y][0]= imageRGBA[idxPixelW][idxPixelH][0];
    }
  }

  // Iteratively smooth the field
  for (int k= 0; k < D.UI[SmoothIter______].I(); k++) {
    std::vector<std::vector<std::vector<double>>> fieldOld= field;
    for (int x= 1; x < nX - 1; x++) {
      for (int y= 1; y < nY - 1; y++) {
        field[x][y][0]= (fieldOld[x][y][0] + fieldOld[x - 1][y][0] + fieldOld[x + 1][y][0] + fieldOld[x][y - 1][0] + fieldOld[x][y + 1][0]) / 5.0f;
      }
    }
  }

  // Replicate the field along extrusion direction and add top/bottom empty layers
  std::vector<std::vector<std::vector<double>>> fieldOld= field;
  field= Field::AllocField3D(nX, nY, 4, 0.0);
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      field[x][y][1]= fieldOld[x][y][0];
      field[x][y][2]= fieldOld[x][y][0];
    }
  }

  // Compute the isosurface with marching cubes
  std::vector<std::array<double, 3>> oVertices;
  std::vector<std::array<int, 3>> oTriangles;
  MarchingCubes::ComputeMarchingCubes(D.UI[Isovalue________].F(), D.boxMin, D.boxMax, field, oVertices, oTriangles);

  // Snap the mesh nodes to the box in the extrusion direction
  for (unsigned int k= 0; k < oVertices.size(); k++) {
    oVertices[k][2]= (oVertices[k][2] < 0.5 * (D.boxMin[2] + D.boxMax[2])) ? D.boxMin[2] : D.boxMax[2];
  }

  // Write the obj file on disk
  std::vector<std::array<double, 3>> oVerticesColors;
  std::vector<std::array<int, 4>> oQuads;
  FileOutput::SaveMeshOBJFile("FileOutput/test.obj", oVertices, oVerticesColors, oTriangles, oQuads, D.UI[VerboseLevel____].I());
}


// Handle keypress
void ImageExtruMesh::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
  if (D.UI[VerboseLevel____].I() >= 5) printf("ImageExtruMesh::KeyPress()\n");
}


// Animate the project
void ImageExtruMesh::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
  if (D.UI[VerboseLevel____].I() >= 5) printf("ImageExtruMesh::Animate()\n");
}


// Draw the project
void ImageExtruMesh::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;
  if (D.UI[VerboseLevel____].I() >= 5) printf("ImageExtruMesh::Draw()\n");
}
