#include "ImageExtruMesh.hpp"


// Standard lib
#include <cmath>
#include <cstring>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "FileIO/FileOutput.hpp"
#include "Geom/BoxGrid.hpp"
#include "Geom/MarchingCubes.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"

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
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ResolutionX_____", 100));
    D.UI.push_back(ParamUI("ResolutionY_____", 100));
    D.UI.push_back(ParamUI("ResolutionZ_____", 20));
    D.UI.push_back(ParamUI("SizeX___________", 1.0));
    D.UI.push_back(ParamUI("SizeY___________", 1.0));
    D.UI.push_back(ParamUI("SizeZ___________", 0.2));
    D.UI.push_back(ParamUI("CenterX_________", 0.5));
    D.UI.push_back(ParamUI("CenterY_________", 0.5));
    D.UI.push_back(ParamUI("CenterZ_________", 0.5));
    D.UI.push_back(ParamUI("HeightMapMode___", 1));
    D.UI.push_back(ParamUI("HeightMapSlope__", 0.1));
    D.UI.push_back(ParamUI("ValOffset_______", 0.5));
    D.UI.push_back(ParamUI("BaseRelHeight___", 0.5));
    D.UI.push_back(ParamUI("GeomSnapMode____", 1));
    D.UI.push_back(ParamUI("MidRelHeight____", 0.3));
    D.UI.push_back(ParamUI("SphereShift_____", -10.0));
    D.UI.push_back(ParamUI("OffsetScaling___", 1.0));
    D.UI.push_back(ParamUI("SmoothIter______", 4));
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
  for (int idxParam= ResolutionX_____; idxParam <= Isovalue________; idxParam++)
    if (D.UI[idxParam].hasChanged()) isRefreshed= false;
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
  D.boxMin= {D.UI[CenterX_________].D() - 0.5 * D.UI[SizeX___________].D(),
             D.UI[CenterY_________].D() - 0.5 * D.UI[SizeY___________].D(),
             D.UI[CenterZ_________].D() - 0.5 * D.UI[SizeZ___________].D()};
  D.boxMax= {D.UI[CenterX_________].D() + 0.5 * D.UI[SizeX___________].D(),
             D.UI[CenterY_________].D() + 0.5 * D.UI[SizeY___________].D(),
             D.UI[CenterZ_________].D() + 0.5 * D.UI[SizeZ___________].D()};

  // Get dimensions
  const int nX= std::max(D.UI[ResolutionX_____].I(), 1);
  const int nY= std::max(D.UI[ResolutionY_____].I(), 1);
  const int nZ= std::max(D.UI[ResolutionZ_____].I(), 1);

  // Load bitmap file and get dimensions
  int nW, nH;
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  FileInput::LoadImageBMPFile("./FileInput/Images/Logo.bmp", imageRGBA, false);
  Field::GetFieldDimensions(imageRGBA, nW, nH);
  if (nW < 1 || nH < 1) return;

  // Sample field from image with zero order interpolation (nearest neighbor) scheme
  ScalarField= Field::AllocField3D(nX, nY, nZ, 0.0);
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      const float posW= (float)(nW - 1) * ((float)x + 0.5f) / (float)nX;
      const float posH= (float)(nH - 1) * ((float)y + 0.5f) / (float)nY;
      const int idxPixelW= std::min(std::max((int)std::round(posW), 0), nW - 1);
      const int idxPixelH= std::min(std::max((int)std::round(posH), 0), nH - 1);
      ScalarField[x][y][0]= (imageRGBA[idxPixelW][idxPixelH][0] + imageRGBA[idxPixelW][idxPixelH][1] + imageRGBA[idxPixelW][idxPixelH][2]) / 3.0;
    }
  }

  // Spread value along Z axis of the field
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 1; z < nZ; z++) {
        ScalarField[x][y][z]= ScalarField[x][y][z - 1] - std::max(D.UI[HeightMapSlope__].D(), 0.0);
      }
    }
  }

  // Apply the value offset
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 1; z < nZ; z++)
        ScalarField[x][y][z]+= D.UI[ValOffset_______].D();

  // Add the optional flat circular base
  if (D.UI[BaseRelHeight___].D() > 0.0 && nZ > 1) {
    for (int x= 1; x < nX - 1; x++)
      for (int y= 1; y < nY - 1; y++)
        for (int z= 1; z < nZ - 1; z++)
          if (std::pow(double(x) / double(nX - 1) - 0.5, 2.0) + std::pow(double(y) / double(nY - 1) - 0.5, 2.0) < 0.25)
            if (double(z) / double(nZ - 1) < D.UI[BaseRelHeight___].D())
              ScalarField[x][y][z]= 1.0;
  }

  // Iteratively smooth the field
  for (int k= 0; k < D.UI[SmoothIter______].I(); k++) {
    std::vector<std::vector<std::vector<double>>> fieldOld= ScalarField;
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          ScalarField[x][y][z]= fieldOld[x][y][z];
          if (x - 1 >= 0) ScalarField[x][y][z]+= fieldOld[x - 1][y][z];
          if (x + 1 < nX) ScalarField[x][y][z]+= fieldOld[x + 1][y][z];
          if (y - 1 >= 0) ScalarField[x][y][z]+= fieldOld[x][y - 1][z];
          if (y + 1 < nY) ScalarField[x][y][z]+= fieldOld[x][y + 1][z];
          if (z - 1 >= 0) ScalarField[x][y][z]+= fieldOld[x][y][z - 1];
          if (z + 1 < nZ) ScalarField[x][y][z]+= fieldOld[x][y][z + 1];
          ScalarField[x][y][z]/= 7.0;
        }
      }
    }
  }

  // Compute the isosurface with marching cubes on expanded field to avoid border holes
  double stepX, stepY, stepZ;
  BoxGrid::GetVoxelSizes(nX, nY, nZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ);
  std::array<double, 3> tmpBoxMin= {D.boxMin[0] - stepX, D.boxMin[1] - stepY, D.boxMin[2] - stepZ};
  std::array<double, 3> tmpBoxMax= {D.boxMax[0] + stepX, D.boxMax[1] + stepY, D.boxMax[2] + stepZ};
  std::vector<std::vector<std::vector<double>>> TmpField= Field::AllocField3D(nX + 2, nY + 2, nZ + 2, 0.0);
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        TmpField[x + 1][y + 1][z + 1]= ScalarField[x][y][z];
  Verts.clear();
  VertsCol.clear();
  Tris.clear();
  Quads.clear();
  MarchingCubes::ComputeMarchingCubes(D.UI[Isovalue________].F(), tmpBoxMin, tmpBoxMax, TmpField, Verts, Tris);

  if (D.UI[GeomSnapMode____].I() == 1) {
    // Snap the mesh nodes to the box in the extrusion direction
    for (unsigned int k= 0; k < Verts.size(); k++) {
      Verts[k][2]= (Verts[k][2] < D.boxMin[2] + 0.5 * (D.boxMax[2] - D.boxMin[2])) ? D.boxMin[2] : D.boxMax[2];
    }
  }
  else if (D.UI[GeomSnapMode____].I() == 2) {
    // Snap the mesh nodes to a sphere
    for (unsigned int k= 0; k < Verts.size(); k++) {
      // Snap the bottom mesh nodes to the box in the extrusion direction
      if (Verts[k][2] < D.boxMin[2] + D.UI[MidRelHeight____].D() * (D.boxMax[2] - D.boxMin[2])) {
        Verts[k][2]= D.boxMin[2];
      }
      else {
        // Snap the top mesh nodes to the surface of a sphere
        Vec::Vec3<double> vert(Verts[k][0], Verts[k][1], Verts[k][2]);
        Vec::Vec3<double> sphere(0.5 * (D.boxMin[0] + D.boxMax[0]), 0.5 * (D.boxMin[1] + D.boxMax[1]), D.boxMin[2] + D.UI[SphereShift_____].D());
        double radius= (sphere - Vec::Vec3<double>(D.boxMin[0], D.boxMin[1], D.boxMin[2])).norm();
        double offset= vert[2] - D.UI[MidRelHeight____].D() * (D.boxMax[2] - D.boxMin[2]);
        vert= sphere + (vert - sphere).normalized() * (radius + offset * D.UI[OffsetScaling___].D());
        Verts[k]= std::array<double, 3>{vert[0], vert[1], vert[2]};
      }
    }
  }
}


// Handle keypress
void ImageExtruMesh::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
  if (D.UI[VerboseLevel____].I() >= 5) printf("ImageExtruMesh::KeyPress()\n");

  // Write the obj file on disk
  FileOutput::SaveMeshOBJFile("FileOutput/test.obj", Verts, VertsCol, Tris, Quads, D.UI[VerboseLevel____].I());
}


// Handle mouse action
void ImageExtruMesh::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
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

  // Draw triangles
  if (D.displayMode3) {
    glEnable(GL_LIGHTING);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_TRIANGLES);
    for (int k= 0; k < (int)Tris.size(); k++) {
      Vec::Vec3<double> v0(Verts[Tris[k][0]][0], Verts[Tris[k][0]][1], Verts[Tris[k][0]][2]);
      Vec::Vec3<double> v1(Verts[Tris[k][1]][0], Verts[Tris[k][1]][1], Verts[Tris[k][1]][2]);
      Vec::Vec3<double> v2(Verts[Tris[k][2]][0], Verts[Tris[k][2]][1], Verts[Tris[k][2]][2]);
      Vec::Vec3<double> n0= (v1 - v0).cross(v2 - v0).normalized();
      Vec::Vec3<double> n1= n0;
      Vec::Vec3<double> n2= n0;
      glNormal3f(n0[0], n0[1], n0[2]);
      glVertex3f(v0[0], v0[1], v0[2]);
      glNormal3f(n1[0], n1[1], n1[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
      glNormal3f(n2[0], n2[1], n2[2]);
      glVertex3f(v2[0], v2[1], v2[2]);
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }


  // Draw scalar field
  if (D.displayMode4) {
    int nX, nY, nZ;
    Field::GetFieldDimensions(ScalarField, nX, nY, nZ);
    if (nX > 0 && nY > 0 && nZ > 0) {
      double stepX, stepY, stepZ, voxDiag, startX, startY, startZ;
      BoxGrid::GetVoxelSizes(nX, nY, nZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ, voxDiag);
      BoxGrid::GetVoxelStart(D.boxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      for (int x= 0; x < nX; x++) {
        for (int y= 0; y < nY; y++) {
          for (int z= 0; z < nZ; z++) {
            if (ScalarField[x][y][z] < D.UI[Isovalue________].D()) continue;
            float r= 0.0f, g= 0.0f, b= 0.0f;
            Colormap::RatioToJetBrightSmooth(0.5 * (ScalarField[x][y][z] + 0.5), r, g, b);
            glColor3f(r, g, b);
            glVertex3f(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
          }
        }
      }
      glEnd();
      glPointSize(1.0f);
    }
  }
}
