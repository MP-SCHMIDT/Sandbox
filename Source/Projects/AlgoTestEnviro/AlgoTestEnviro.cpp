#include "AlgoTestEnviro.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <limits>
#include <string>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "FileIO/FileOutput.hpp"
#include "Geom/BoxGrid.hpp"
#include "Geom/MarchingCubes.hpp"
#include "Geom/MeshVoxelizer.hpp"
#include "Geom/DistanceTransform.hpp"
#include "Geom/Labelizer.hpp"
#include "Geom/PrimitiveCSG.hpp"
#include "Geom/Sketch.hpp"
#include "Type/Field.hpp"
#include "Math/Functions.hpp"
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
AlgoTestEnviro::AlgoTestEnviro() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void AlgoTestEnviro::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("TestParamALG_00_", 0.9));
    D.UI.push_back(ParamUI("TestParamALG_01_", 0.1));
    D.UI.push_back(ParamUI("TestParamALG_02_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_03_", 64.0));
    D.UI.push_back(ParamUI("TestParamALG_04_", 64.0));
    D.UI.push_back(ParamUI("TestParamALG_05_", 64.0));
    D.UI.push_back(ParamUI("TestParamALG_06_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_07_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_08_", 1.e6));
    D.UI.push_back(ParamUI("TestParamALG_09_", 2.0));
    D.UI.push_back(ParamUI("TestParamALG_10_", 2.0));
    D.UI.push_back(ParamUI("TestParamALG_11_", 2.0));
    D.UI.push_back(ParamUI("TestParamALG_12_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_13_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_14_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_15_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_16_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_17_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_18_", 0.0));
    D.UI.push_back(ParamUI("TestParamALG_19_", 0.0));
    D.UI.push_back(ParamUI("Isocut__________", 0.0));
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));
    D.UI.push_back(ParamUI("VerboseLevel____", 1));

    D.displayModeLabel[1]= "Vert";
    D.displayModeLabel[2]= "Bar";
    D.displayModeLabel[3]= "Tri";
    D.displayModeLabel[4]= "Field";
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool AlgoTestEnviro::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool AlgoTestEnviro::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void AlgoTestEnviro::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Allocate()\n");
}


// Refresh the project
void AlgoTestEnviro::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  if (D.UI[VerboseLevel____].I() >= 5) printf("Refresh()\n");
}


// Handle UI parameter change
void AlgoTestEnviro::ParamChange() {
  if (D.UI[VerboseLevel____].I() >= 5) printf("ParamChange()\n");
}


// Handle keypress
void AlgoTestEnviro::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();

  if (D.UI[VerboseLevel____].I() >= 5) printf("KeyPress()\n");

  // Testing CSG field
  if (D.keyLetterUpperCase == 'F') {
    int const nX= D.UI[TestParamALG_03_].I();
    int const nY= D.UI[TestParamALG_04_].I();
    int const nZ= D.UI[TestParamALG_05_].I();
    ScalarField= Field::Field3(nX, nY, nZ, std::numeric_limits<double>::max());
    std::array<double, 3> P0({0.1, 0.2, 0.3});
    std::array<double, 3> P1({0.8, 0.6, 0.7});
    PrimitiveCSG::ConeRound(nX, nY, nZ, P0, P1, 0.1, 0.2, PrimitiveCSG::BooleanMode::Union, D.boxMin, D.boxMax, ScalarField.data);
  }

  // Testing sketch smoothing
  if (D.keyLetterUpperCase == 'S') {
    std::vector<std::array<double, 3>> polylineRef;
    polylineRef.push_back(std::array<double, 3>({0.0, 0.0, 0.0}));
    polylineRef.push_back(std::array<double, 3>({0.5, 0.5, 0.5}));
    polylineRef.push_back(std::array<double, 3>({1.0, -0.5, -0.5}));
    std::vector<std::array<double, 3>> polylineNew= polylineRef;
    Sketch::PolylineSubdivideAndSmooth(false, 3, 2, polylineNew);

    Verts.clear();
    Bars.clear();
    Tris.clear();
    for (int k= 0; k < (int)polylineNew.size(); k++) {
      Verts.push_back(polylineNew[k]);
      Bars.push_back(std::array<int, 2>{k, (k + 1) % (int)polylineNew.size()});
    }
  }

  // Input OBJ mesh
  if (D.keyLetterUpperCase == 'L') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    Verts.clear();
    Bars.clear();
    Tris.clear();
    std::vector<std::array<double, 3>> VertsCol;
    FileInput::LoadMeshOBJFile("FileOutput/Test.obj", Verts, VertsCol, Tris, D.UI[VerboseLevel____].I());
    if (D.UI[VerboseLevel____].I() >= 1) printf("LoadMesh %.3f ", Timer::PopTimer());
  }

  // Output OBJ mesh
  if (D.keyLetterUpperCase == 'O') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    std::vector<std::array<double, 3>> VertsCol;
    std::vector<std::array<int, 4>> Quads;
    FileOutput::SaveMeshOBJFile("FileOutput/TestOut.obj", Verts, VertsCol, Tris, Quads, D.UI[VerboseLevel____].I());
    if (D.UI[VerboseLevel____].I() >= 1) printf("SaveMesh %.3f ", Timer::PopTimer());
  }

  // Marching Cubes on current scalar field
  if (D.keyLetterUpperCase == 'M') {
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    Verts.clear();
    Bars.clear();
    Tris.clear();
    MarchingCubes::BuildMesh(ScalarField.nX, ScalarField.nY, ScalarField.nZ, 0.0, true, false, D.boxMin, D.boxMax, ScalarField.data, Verts, Tris);
    if (D.UI[VerboseLevel____].I() >= 1) printf("MarchingCubes %.3f ", Timer::PopTimer());
  }

  // Build clean SDF from current unclean mesh
  if (D.keyLetterUpperCase == 'R') {
    if (Verts.size() <= 2) return;
    if (Tris.size() <= 0) return;

    const int nbVox= D.UI[TestParamALG_08_].I();           // Desired number of voxels
    const int nbMargin= D.UI[TestParamALG_09_].I();        // Number of voxels making the margin around the model
    const float nbVoxDilation= D.UI[TestParamALG_10_].F(); // Distance offset to dilate the SDF and close holes
    const float nbVoxErosion= D.UI[TestParamALG_11_].F();  // Distance offset to erode the SDF and move the isolevel back to the original location

    // Compute the mesh bounding box
    D.boxMin= {Verts[0][0], Verts[0][1], Verts[0][2]};
    D.boxMax= {Verts[0][0], Verts[0][1], Verts[0][2]};
    for (unsigned int k= 1; k < Verts.size(); k++) {
      for (unsigned int dim= 0; dim < 3; dim++) {
        if (D.boxMin[dim] > Verts[k][dim]) D.boxMin[dim]= Verts[k][dim];
        if (D.boxMax[dim] < Verts[k][dim]) D.boxMax[dim]= Verts[k][dim];
      }
    }
  
    // Compute the grid resolution with margin
    const double boxDX= D.boxMax[0] - D.boxMin[0];
    const double boxDY= D.boxMax[1] - D.boxMin[1];
    const double boxDZ= D.boxMax[2] - D.boxMin[2];
    const double stepSize= std::pow(boxDX * boxDY * boxDZ, 1.0/3.0) / std::pow((double)nbVox, 1.0/3.0);
    const int nX= std::ceil(boxDX / stepSize) + 2 * nbMargin;
    const int nY= std::ceil(boxDY / stepSize) + 2 * nbMargin;
    const int nZ= std::ceil(boxDZ / stepSize) + 2 * nbMargin;
    const int nXYZ= nX*nY*nZ;

    std::array<double, 3> boxMinOld= D.boxMin;
    std::array<double, 3> boxMaxOld= D.boxMax;
    D.boxMin[0]= 0.5*(boxMinOld[0]+boxMaxOld[0]) - 0.5*double(nX+1)*stepSize;
    D.boxMin[1]= 0.5*(boxMinOld[1]+boxMaxOld[1]) - 0.5*double(nY+1)*stepSize;
    D.boxMin[2]= 0.5*(boxMinOld[2]+boxMaxOld[2]) - 0.5*double(nZ+1)*stepSize;
    D.boxMax[0]= 0.5*(boxMinOld[0]+boxMaxOld[0]) + 0.5*double(nX+1)*stepSize;
    D.boxMax[1]= 0.5*(boxMinOld[1]+boxMaxOld[1]) + 0.5*double(nY+1)*stepSize;
    D.boxMax[2]= 0.5*(boxMinOld[2]+boxMaxOld[2]) + 0.5*double(nZ+1)*stepSize;
  
    ScalarField= Field::Field3<double>(nX, nY, nZ, 0.0);

    // Compute the narrow band signed distance field
    Field::Field3<char> isFixed;
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    MeshVoxelizer::ComputeNarrowBandSignedDistanceField(Verts, Tris, D.boxMin, D.boxMax, nX, nY, nZ, true, isFixed, ScalarField);
    if (D.UI[VerboseLevel____].I() >= 1) printf("NBSDF %.3f ", Timer::PopTimer());

    // Spread signed distances to the full domain
    Field::Field3<float> SolidSDF(nX, nY, nZ, 0.0f);
    for (int xyz= 0; xyz < nXYZ; xyz++)
      SolidSDF.at(xyz)= ScalarField.at(xyz);
    if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
    DistanceTransform::ApplyDistanceTransform(1, stepSize, stepSize, stepSize, isFixed, SolidSDF, 0);
    if (D.UI[VerboseLevel____].I() >= 1) printf("DistTransf %.3f ", Timer::PopTimer());
  
    // Repair unreliable signs of SDF on non solid surface meshes to close holes and delete internal geometries
    if (nbVoxDilation > 0.0f) {
      // Apply absolute value and dilation offset
      for (int xyz= 0; xyz < nXYZ; xyz++)
        SolidSDF.at(xyz)= std::abs(SolidSDF.at(xyz)) - nbVoxDilation * stepSize;
  
      // Label the regions based on distance sign
      if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
      Field::Field3<int> category(nX, nY, nZ, 0);
      Field::Field3<int> component(nX, nY, nZ, 0);
      for (int xyz= 0; xyz < nXYZ; xyz++)
        if (SolidSDF.at(xyz) < 0.0f)
          category.at(xyz)= 1;
      Labelizer::LabelConnectedComponents(nX, nY, nZ, category.data, component.data);
      if (D.UI[VerboseLevel____].I() >= 1) printf("Label %.3f ", Timer::PopTimer());
  
      // Find the main outside label by checking the comain corners
      int refLabel= component.at(1, 1, 1);
      if      (SolidSDF.at(   1,    1, nZ-1) > 0.0f) refLabel= component.at(   1,    1, nZ-1);
      else if (SolidSDF.at(   1, nY-1,    1) > 0.0f) refLabel= component.at(   1, nY-1,    1);
      else if (SolidSDF.at(   1, nY-1, nZ-1) > 0.0f) refLabel= component.at(   1, nY-1, nZ-1);
      else if (SolidSDF.at(nX-1,    1,    1) > 0.0f) refLabel= component.at(nX-1,    1,    1);
      else if (SolidSDF.at(nX-1,    1, nZ-1) > 0.0f) refLabel= component.at(nX-1,    1, nZ-1);
      else if (SolidSDF.at(nX-1, nY-1,    1) > 0.0f) refLabel= component.at(nX-1, nY-1,    1);
      else if (SolidSDF.at(nX-1, nY-1, nZ-1) > 0.0f) refLabel= component.at(nX-1, nY-1, nZ-1);
      
      // Flip the sign of any fluid region not connected to the main outside region
      for (int xyz= 0; xyz < nXYZ; xyz++)
        if (SolidSDF.at(xyz) > 0.0f && component.at(xyz) != refLabel)
          SolidSDF.at(xyz)= -SolidSDF.at(xyz);
  
      // Repair SDF using the offset distance value as reference and re spreading inward before counter-offsetting
      if (nbVoxErosion > 0.0f) {
        // Only keep distances straddling the offset interface
        if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
        isFixed= Field::Field3<char> (nX, nY, nZ, 0);
        for (int x= 0; x < nX; x++) {
          for (int y= 0; y < nY; y++) {
            for (int z= 0; z < nZ; z++) {
              const int xyz= SolidSDF.getFlatIndex(x, y, z);
              for (int xOff= std::max(x - 1, 0); xOff <= std::min(x + 1, nX - 1); xOff++) {
                for (int yOff= std::max(y - 1, 0); yOff <= std::min(y + 1, nY - 1); yOff++) {
                  for (int zOff= std::max(z - 1, 0); zOff <= std::min(z + 1, nZ - 1); zOff++) {
                    const int xyzOff= SolidSDF.getFlatIndex(xOff, yOff, zOff);
                    if (SolidSDF.at(xyz) * SolidSDF.at(xyzOff) <= 0.0f) {
                      isFixed.at(xyz)= 1;
                      isFixed.at(xyzOff)= 1;
                    }
                  }
                }
              }
            }
          }
        }
        if (D.UI[VerboseLevel____].I() >= 1) printf("FlagInterf %.3f ", Timer::PopTimer());
        
        // Spread signed distances to the full domain
        if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();
        DistanceTransform::ApplyDistanceTransform(1, stepSize, stepSize, stepSize, isFixed, SolidSDF, 0);
        if (D.UI[VerboseLevel____].I() >= 1) printf("DistTransf %.3f ", Timer::PopTimer());
  
        // Apply the erosion offset
        for (int xyz= 0; xyz < nXYZ; xyz++)
          SolidSDF.at(xyz)+= nbVoxErosion * stepSize;
      }
    }
    for (int xyz= 0; xyz < nXYZ; xyz++)
      ScalarField.at(xyz)= SolidSDF.at(xyz);
  }

  // Testing Penal function
  if (D.keyLetterUpperCase == 'P') {
    Verts.clear();
    for (int k= 0; k < 1000; k++) {
      const double ratio= double(k) / double(1000 - 1);
      Verts.push_back(std::array<double, 3>{D.boxMin[0] + (D.boxMax[0] - D.boxMin[0]) * 0.5,
                                            D.boxMin[1] + (D.boxMax[1] - D.boxMin[1]) * ratio,
                                            D.boxMin[2] + (D.boxMax[2] - D.boxMin[2]) * Functions::PenalInterpo(ratio, D.UI[TestParamALG_09_].D(), false)});
      Verts.push_back(std::array<double, 3>{D.boxMin[0] + (D.boxMax[0] - D.boxMin[0]) * 0.5,
                                            D.boxMin[1] + (D.boxMax[1] - D.boxMin[1]) * ratio,
                                            D.boxMin[2] + (D.boxMax[2] - D.boxMin[2]) * Functions::PenalInterpo(ratio, D.UI[TestParamALG_09_].D(), true)});
    }
  }
}


// Handle mouse action
void AlgoTestEnviro::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (D.UI[VerboseLevel____].I() >= 5) printf("MousePress()\n");
}


// Animate the project
void AlgoTestEnviro::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();
  if (D.UI[VerboseLevel____].I() >= 5) printf("Animate()\n");
}


// Draw the project
void AlgoTestEnviro::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;
  if (D.UI[VerboseLevel____].I() >= 5) printf("Draw()\n");

  // Draw vertices
  if (D.displayMode[1]) {
    glPointSize(5.0f);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_POINTS);
    for (int k= 0; k < (int)Verts.size(); k++) {
      glVertex3dv(Verts[k].data());
    }
    glEnd();
    glPointSize(1.0f);
  }

  // Draw bars
  if (D.displayMode[2]) {
    glColor3f(0.6f, 0.6f, 0.6f);
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    for (int k= 0; k < (int)Bars.size(); k++) {
      float r= 0.5 + 0.5 * (Verts[Bars[k][0]][0]-Verts[Bars[k][1]][0]);
      float g= 0.5 + 0.5 * (Verts[Bars[k][0]][1]-Verts[Bars[k][1]][1]);
      float b= 0.5 + 0.5 * (Verts[Bars[k][0]][2]-Verts[Bars[k][1]][2]);
      glColor3f(r, g, b);
      glVertex3dv(Verts[Bars[k][0]].data());
      glVertex3dv(Verts[Bars[k][1]].data());
    }
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw triangles
  if (D.displayMode[3]) {
    glEnable(GL_LIGHTING);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_TRIANGLES);
    for (int k= 0; k < (int)Tris.size(); k++) {
      Vec::Vec3<double> v0(Verts[Tris[k][0]]);
      Vec::Vec3<double> v1(Verts[Tris[k][1]]);
      Vec::Vec3<double> v2(Verts[Tris[k][2]]);
      Vec::Vec3<double> n0= (v1 - v0).cross(v2 - v0).normalized();
      Vec::Vec3<double> n1= n0;
      Vec::Vec3<double> n2= n0;
      glNormal3dv(n0.array());
      glVertex3dv(v0.array());
      glNormal3dv(n1.array());
      glVertex3dv(v1.array());
      glNormal3dv(n2.array());
      glVertex3dv(v2.array());
    }
    glEnd();
    glDisable(GL_LIGHTING);
  }

  // Draw scalar field
  if (D.displayMode[4]) {
    if (ScalarField.nXYZ > 0) {
      double stepX, stepY, stepZ, voxDiag, startX, startY, startZ;
      BoxGrid::GetVoxelSizes(ScalarField.nX, ScalarField.nY, ScalarField.nZ, D.boxMin, D.boxMax, true, stepX, stepY, stepZ, voxDiag);
      BoxGrid::GetVoxelStart(D.boxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      for (int x= 0; x < ScalarField.nX; x++) {
        for (int y= 0; y < ScalarField.nY; y++) {
          for (int z= 0; z < ScalarField.nZ; z++) {
            if (ScalarField.at(x, y, z) < D.UI[Isocut__________].D()) continue;
            float r= 0.0f, g= 0.0f, b= 0.0f;
            Colormap::RatioToJetBrightSmooth(0.5 + 0.5 * ScalarField.at(x, y, z) * D.UI[ColorFactor_____].D(), r, g, b);
            glColor3f(r, g, b);
            glVertex3f(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
          }
        }
      }
      glEnd();
      glPointSize(1.0f);
    }
  }

  // Draw the status
  D.Status.clear();
  D.Status.resize(2);
  D.Status[0]= std::format("Grid={} x {} x {} = {}", ScalarField.nX, ScalarField.nY, ScalarField.nZ, ScalarField.nXYZ);
  D.Status[1]= std::format("Verts={} Bars={} Tris={}", Verts.size(), Bars.size(), Tris.size());
}
