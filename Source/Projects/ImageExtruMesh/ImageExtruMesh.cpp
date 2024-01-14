#include "ImageExtruMesh.hpp"


// Standard lib
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "FileIO/FileInput.hpp"
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
    D.UI.push_back(ParamUI("InputFile_______", 1));
    D.UI.push_back(ParamUI("DomainW_________", 50));
    D.UI.push_back(ParamUI("DomainH_________", 50));
    D.UI.push_back(ParamUI("SmoothIter______", 3));
    D.UI.push_back(ParamUI("Isovalue________", 0.5));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));
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

  // Get and check dimensions
  int nX= 1;
  int nY= D.UI[DomainW_________].I();
  int nZ= D.UI[DomainH_________].I();
  if (nX < 1 || nY < 1 || nZ < 1) return;

  // Load bitmap file
  const int inputFile= D.UI[InputFile_______].I();
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (inputFile == 0) FileInput::LoadImageBMPFile("./FileInput/Images/Logo.bmp", imageRGBA, false);
  if (imageRGBA.empty()) return;

  // Project bitmap alpha channel as scalar field values
  std::vector<std::vector<std::vector<double>>> field= Field::AllocField3D(1, nY, nZ, 0.0);
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        const float posW= (float)(imageRGBA.size() - 1) * ((float)y + 0.5f) / (float)nY;
        const float posH= (float)(imageRGBA[0].size() - 1) * ((float)z + 0.5f) / (float)nZ;
        const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
        const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
        const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
        field[x][y][z]= 1.0 - colRGBA[3];
      }
    }
  }

  // Iteratively smooth the field
  for (int k= 0; k < D.UI[SmoothIter______].I(); k++) {
    std::vector<std::vector<std::vector<double>>> fieldOld= field;
    for (int y= 1; y < nY - 1; y++) {
      for (int z= 1; z < nZ - 1; z++) {
        field[0][y][z]= (fieldOld[0][y][z] + fieldOld[0][y + 1][z] + fieldOld[0][y - 1][z] + fieldOld[0][y][z + 1] + fieldOld[0][y][z - 1]) / 5.0f;
      }
    }
  }

  // Replicate the field along extrusion direction and add top/bottom empty layers
  std::vector<std::vector<std::vector<double>>> fieldOld= field;
  field= Field::AllocField3D(nX + 4, nY, nZ, 0.0);
  field[1]= fieldOld[0];
  field[2]= fieldOld[0];

  // Compute the isosurface with marching cubes
  std::vector<std::array<double, 3>> oVertices;
  std::vector<std::array<int, 3>> oTriangles;
  MarchingCubes::ComputeMarchingCubes(0.3f, std::array<double, 3>({0.0, 0.0, 0.0}), std::array<double, 3>({1.0, 1.0, 1.0}), field, oVertices, oTriangles);

  // Snap the mesh nodes to the box in the extrusion direction
  for (unsigned int k= 0; k < oVertices.size(); k++) {
    oVertices[k][0]= (oVertices[k][0] < 0.5) ? 0.0 : 1.0;
  }

  // Write the obj file on disk
  std::string iFullpath= "FileOutput/test.obj";
  printf("Saving OBJ mesh file [%s]\n", iFullpath.c_str());
  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n\n");
  }
  else {
    for (unsigned int k= 0; k < oVertices.size(); k++) {
      fprintf(outputFile, "v %lf %lf %lf\n", oVertices[k][0], oVertices[k][1], oVertices[k][2]);
    }
    for (unsigned int k= 0; k < oTriangles.size(); k++) {
      fprintf(outputFile, "f %d %d %d\n", oTriangles[k][0] + 1, oTriangles[k][1] + 1, oTriangles[k][2] + 1);
    }
    fclose(outputFile);
  }
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
