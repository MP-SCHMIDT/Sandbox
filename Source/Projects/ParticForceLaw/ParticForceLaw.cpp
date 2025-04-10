#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <format>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
ParticForceLaw::ParticForceLaw() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void ParticForceLaw::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("DomainX_________", 1.0));    // Simulation domain size
    D.UI.push_back(ParamUI("DomainY_________", 1.0));    // Simulation domain size
    D.UI.push_back(ParamUI("DomainZ_________", 1.0));    // Simulation domain size
    D.UI.push_back(ParamUI("ScenarioPreset__", 4));      // Scenario choice, 0 load file, 1+ hard coded scenarios
    D.UI.push_back(ParamUI("Scenario2DID____", 0));      // BMP file to load
    D.UI.push_back(ParamUI("Scenario2DThick_", 0.5));    // Relative thickness of loaded file wrt domain
    D.UI.push_back(ParamUI("LatticePitch____", 0.04));   // Charateristic lengthscale of particle cloud, controls the number of particles
    D.UI.push_back(ParamUI("LatticePattern__", 2));      // Choose the particle cloud pattern
    D.UI.push_back(ParamUI("ConstrainDim2D__", 0));      // Optionally constrain movement to plane orthogonal to X=1, Y=2, X=3 axis
    D.UI.push_back(ParamUI("______________00", NAN));    //
    D.UI.push_back(ParamUI("StepsPerDraw____", 1));      // Number of simulation steps between each draw call
    D.UI.push_back(ParamUI("TimeStep________", 1.e-3));  // Physical timestep of explicit integration
    D.UI.push_back(ParamUI("BucketCapacity__", 30));     // Chosen bucket capacity
    D.UI.push_back(ParamUI("BucketMaxCount__", 100000)); // Maximum number of buckets allocated for the spatial partition
    D.UI.push_back(ParamUI("BucketFitDomain_", 1));      // Flag to fit the buckets to the particle cloud bounding box
    D.UI.push_back(ParamUI("DampingRadRel___", 0.1));    // Physical interparticle damping proportional to difference of radial velocities
    D.UI.push_back(ParamUI("DampingVelRel___", 0.0));    // Unphysical absolute velocity damping for stability
    D.UI.push_back(ParamUI("MaterialDensity_", 200.0));  // Density of simulated material
    D.UI.push_back(ParamUI("______________01", NAN));    //
    D.UI.push_back(ParamUI("BCVelX__________", 0.0));    // Velocity for particles with associated BC
    D.UI.push_back(ParamUI("BCVelY__________", 0.0));    // Velocity for particles with associated BC
    D.UI.push_back(ParamUI("BCVelZ__________", 1.0));    // Velocity for particles with associated BC
    D.UI.push_back(ParamUI("BCForX__________", 0.0));    // External force for particles with associated BC
    D.UI.push_back(ParamUI("BCForY__________", 0.0));    // External force for particles with associated BC
    D.UI.push_back(ParamUI("BCForZ__________", 1.0));    // External force for particles with associated BC
    D.UI.push_back(ParamUI("ForceControl____", 2));      // Force controller mode instead of simple overwrite
    D.UI.push_back(ParamUI("BCPosCoeff______", 1.0));    // Scaling coeff for force controller on position BC
    D.UI.push_back(ParamUI("BCVelCoeff______", 1.0));    // Scaling coeff for force controller on velocity BC
    D.UI.push_back(ParamUI("______________02", NAN));    //
    D.UI.push_back(ParamUI("ForceLawPresetA_", 0));      // Chosen force law for material 1
    D.UI.push_back(ParamUI("ForceLawPresetB_", 1));      // Chosen force law for material 2
    D.UI.push_back(ParamUI("ForceLawScaleA__", 100.0));  // Force law scaling for material 1
    D.UI.push_back(ParamUI("ForceLawScaleB__", 100.0));  // Force law scaling for material 2
    D.UI.push_back(ParamUI("______________03", NAN));    //
    D.UI.push_back(ParamUI("ForceLawA_0_00__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_0_80__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_0_90__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_0_95__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_00__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_05__", -1.0));   // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_10__", -1.0));   // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_20__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_30__", 1.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_40__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_1_50__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_2_00__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_2_50__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("ForceLawA_3_00__", 0.0));    // Force law control points for material 1
    D.UI.push_back(ParamUI("______________04", NAN));    //
    D.UI.push_back(ParamUI("ForceLawB_0_00__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_0_80__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_0_90__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_0_95__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_00__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_05__", -1.0));   // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_10__", -1.0));   // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_20__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_30__", 1.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_40__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_1_50__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_2_00__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_2_50__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("ForceLawB_3_00__", 0.0));    // Force law control points for material 2
    D.UI.push_back(ParamUI("______________05", NAN));    //
    D.UI.push_back(ParamUI("MetaballVoxSize_", 0.04));   // Voxel size for Metaball and Marching Cubes isosurface
    D.UI.push_back(ParamUI("MetaballIsoval__", 0.5));    // Isovalue for Marching Cubes
    D.UI.push_back(ParamUI("ColorMode_______", 3));      // Show material, boundary conditions, velocity, force magnitude, velocity magnitude
    D.UI.push_back(ParamUI("ColorFactor_____", 1.0));    // Color factor for modes supporting it
    D.UI.push_back(ParamUI("VisuScale_______", 0.5));    // Size scaling of display elements
    D.UI.push_back(ParamUI("VisuSimple______", 0));      // Toggle for simplified draw mode
    D.UI.push_back(ParamUI("VisuHideOOB_____", 0));      // Hide particles out of bounds
    D.UI.push_back(ParamUI("VisuMinNeighbor_", 0));      // Hide particles with too few neighbors
    D.UI.push_back(ParamUI("______________06", NAN));    //
    D.UI.push_back(ParamUI("TestParamMIP_0__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_1__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_2__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_3__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_4__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("TestParamMIP_5__", 0));      // Generic param for testing purposes
    D.UI.push_back(ParamUI("VerboseLevel____", 0));      // Verbose mode

    D.displayModeLabel[1]= "Partic";
    D.displayModeLabel[2]= "SpringsBC";
    D.displayModeLabel[3]= "Bucket";
    D.displayModeLabel[4]= "Tris";
    D.displayMode[3]= false;
    D.displayMode[4]= false;
  }

  if (D.UI.size() != VerboseLevel____ + 1) printf("[ERROR] Invalid parameter count in UI\n");

  RunID= 0;

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool ParticForceLaw::CheckAlloc() {
  if (D.UI[DomainX_________].hasChanged()) isAllocated= false;
  if (D.UI[DomainY_________].hasChanged()) isAllocated= false;
  if (D.UI[DomainZ_________].hasChanged()) isAllocated= false;
  if (D.UI[ScenarioPreset__].hasChanged()) isAllocated= false;
  if (D.UI[Scenario2DID____].hasChanged()) isAllocated= false;
  if (D.UI[Scenario2DThick_].hasChanged()) isAllocated= false;
  if (D.UI[LatticePitch____].hasChanged()) isAllocated= false;
  if (D.UI[LatticePattern__].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool ParticForceLaw::CheckRefresh() {
  if (D.UI[ForceLawPresetA_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawPresetB_].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawScaleA__].hasChanged()) isRefreshed= false;
  if (D.UI[ForceLawScaleB__].hasChanged()) isRefreshed= false;
  for (int idxParam= ForceLawA_0_00__; idxParam <= ForceLawA_3_00__; idxParam++)
    if (D.UI[idxParam].hasChanged()) isRefreshed= false;
  for (int idxParam= ForceLawB_0_00__; idxParam <= ForceLawB_3_00__; idxParam++)
    if (D.UI[idxParam].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void ParticForceLaw::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Reset data arrays
  Ref.clear();
  Pos.clear();
  Vel.clear();
  For.clear();
  Col.clear();
  Mat.clear();
  ForceMag.clear();
  Neighbors.clear();
  Sensor.clear();
  BCPos.clear();
  BCVel.clear();
  BCFor.clear();
  MetaballIsUpdated= false;
  Verts.clear();
  Tris.clear();

  // Get domain dimensions
  D.boxMin= {0.5 - 0.5 * std::max(D.UI[DomainX_________].D(), D.UI[LatticePitch____].D()),
             0.5 - 0.5 * std::max(D.UI[DomainY_________].D(), D.UI[LatticePitch____].D()),
             0.5 - 0.5 * std::max(D.UI[DomainZ_________].D(), D.UI[LatticePitch____].D())};
  D.boxMax= {0.5 + 0.5 * std::max(D.UI[DomainX_________].D(), D.UI[LatticePitch____].D()),
             0.5 + 0.5 * std::max(D.UI[DomainY_________].D(), D.UI[LatticePitch____].D()),
             0.5 + 0.5 * std::max(D.UI[DomainZ_________].D(), D.UI[LatticePitch____].D())};

  // Generate the full point cloud over the domain
  std::vector<Vec::Vec3<float>> pointCloud;
  BuildBaseCloud(pointCloud);

  // Generate the scenario
  BuildScenario(pointCloud);
}


// Refresh the project
void ParticForceLaw::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Generate the force law
  BuildForceLaws();

  // Plot the force laws
  if (D.Plot.size() < 6) D.Plot.resize(6);
  for (int idxMat= 0; idxMat < (int)ForceLaw.size(); idxMat++) {
    D.Plot[3 + idxMat].name= "ForceLaw";
    D.Plot[3 + idxMat].isSymmetric= true;
    if (idxMat > 0) D.Plot[3 + idxMat].isSameRange= true;
    D.Plot[3 + idxMat].val.resize(ForceLaw[idxMat].size());
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      D.Plot[3 + idxMat].val[k]= ForceLaw[idxMat][k];
  }
}


// Handle UI parameter change
void ParticForceLaw::ParamChange() {
}


// Handle keypress
void ParticForceLaw::KeyPress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Handle mouse action
void ParticForceLaw::MousePress() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
}


// Animate the project
void ParticForceLaw::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Step forward simulation with explicit numerical integration
  for (int stepIdx= 0; stepIdx < D.UI[StepsPerDraw____].I(); stepIdx++) {
    StepSimulation();
    MetaballIsUpdated= false;
  }

  // Plot spatial partition occupancy
  if (D.Plot.size() < 1) D.Plot.resize(1);
  D.Plot[0].name= "Buckets";
  D.Plot[0].val.resize(Bucket3D.nXYZ + 2);
  D.Plot[0].val[Bucket3D.nXYZ + 0]= 0;
  D.Plot[0].val[Bucket3D.nXYZ + 1]= Bucket3D.bucketCapacity;
  for (int xyz= 0; xyz < Bucket3D.nXYZ; xyz++)
    D.Plot[0].val[xyz]= (double)Bucket3D.buckets[xyz].size();

  // Plot the kinetic and potential energies
  if (D.Plot.size() < 3) D.Plot.resize(3);
  D.Plot[1].name= "KE";
  D.Plot[2].name= "PE";
  if (D.Plot[1].val.size() < 5000) {
    D.Plot[1].val.reserve(5000);
    D.Plot[2].val.reserve(5000);
    D.Plot[1].val.push_back(0.0f);
    D.Plot[2].val.push_back(0.0f);
    for (int k= 0; k < (int)Pos.size(); k++) {
      D.Plot[1].val[D.Plot[1].val.size() - 1]+= Vel[k].normSquared();
      D.Plot[2].val[D.Plot[2].val.size() - 1]+= ForceMag[k];
    }
  }
}


// Draw the project
void ParticForceLaw::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Set particle color
  for (int k= 0; k < (int)Pos.size(); k++) {
    float r= 0.5f, g= 0.5f, b= 0.5f;
    if (D.UI[ColorMode_______].I() == 0) {
      if (Mat[k] == 0) r= 1.0f;
      if (Mat[k] == 1) g= 1.0f;
      if (Mat[k] == 2) b= 1.0f;
    }
    if (D.UI[ColorMode_______].I() == 1) {
      r= 0.5f + 0.3f * (float)BCPos[k];
      g= 0.5f + 0.3f * (float)BCVel[k];
      b= 0.5f + 0.3f * (float)BCFor[k];
    }
    if (D.UI[ColorMode_______].I() == 2) {
      r= D.UI[ColorFactor_____].F() * Vel[k][0] + 0.5f;
      g= D.UI[ColorFactor_____].F() * Vel[k][1] + 0.5f;
      b= D.UI[ColorFactor_____].F() * Vel[k][2] + 0.5f;
    }
    if (D.UI[ColorMode_______].I() == 3) {
      Colormap::RatioToJetBrightSmooth(ForceMag[k] * D.UI[ColorFactor_____].F(), r, g, b);
    }
    if (D.UI[ColorMode_______].I() == 4) {
      Colormap::RatioToBlackBody(Vel[k].norm() * D.UI[ColorFactor_____].F(), r, g, b);
    }
    Col[k].set(r, g, b);
  }

  // Display particles
  if (D.displayMode[1]) {
    if (D.UI[VisuSimple______].I() > 0) {
      glPointSize(1000.0f * D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F());
      glBegin(GL_POINTS);
      for (int k= 0; k < (int)Pos.size(); k++) {
        if (D.UI[VisuHideOOB_____].I() > 0 &&
            (Pos[k][0] < D.boxMin[0] || Pos[k][0] > D.boxMax[0] ||
             Pos[k][1] < D.boxMin[1] || Pos[k][1] > D.boxMax[1] ||
             Pos[k][2] < D.boxMin[2] || Pos[k][2] > D.boxMax[2])) continue;
        if (Neighbors[k] < D.UI[VisuMinNeighbor_].I()) continue;
        glColor3fv(Col[k].array());
        glVertex3fv(Pos[k].array());
      }
      glEnd();
    }
    else {
      glEnable(GL_LIGHTING);
      for (int k= 0; k < (int)Pos.size(); k++) {
        if (D.UI[VisuHideOOB_____].I() > 0 &&
            (Pos[k][0] < D.boxMin[0] || Pos[k][0] > D.boxMax[0] ||
             Pos[k][1] < D.boxMin[1] || Pos[k][1] > D.boxMax[1] ||
             Pos[k][2] < D.boxMin[2] || Pos[k][2] > D.boxMax[2])) continue;
        if (Neighbors[k] < D.UI[VisuMinNeighbor_].I()) continue;
        glColor3fv(Col[k].array());
        glPushMatrix();
        glTranslatef(Pos[k][0], Pos[k][1], Pos[k][2]);
        glutSolidSphere(D.UI[LatticePitch____].F() * D.UI[VisuScale_______].F(), 12, 6);
        glPopMatrix();
      }
      glDisable(GL_LIGHTING);
    }
  }

  // Display BC forces
  if (D.displayMode[2]) {
    if (D.UI[ForceControl____].I() > 0) {
      glLineWidth(2.0f);
      glPointSize(4.0f);
      glBegin(GL_LINES);
      glColor3f(0.5f, 0.5f, 0.5f);
      for (int k= 0; k < (int)Ref.size(); k++) {
        if (BCPos[k] != 0 || BCVel[k] != 0) {
          glVertex3fv(Pos[k].array());
          glVertex3fv(Ref[k].array());
        }
      }
      glEnd();
      glBegin(GL_POINTS);
      glColor3f(0.5f, 0.5f, 0.5f);
      for (int k= 0; k < (int)Ref.size(); k++)
        if (BCPos[k] != 0 || BCVel[k] != 0)
          glVertex3fv(Ref[k].array());
      glEnd();
      glLineWidth(1.0f);
    }
  }

  // Display spatial partition buckets status
  if (D.displayMode[3]) {
    glLineWidth(2.0f);
    // Get dimensions
    const float stepX= (Bucket3D.boxMax[0] - Bucket3D.boxMin[0]) / float(Bucket3D.nX);
    const float stepY= (Bucket3D.boxMax[1] - Bucket3D.boxMin[1]) / float(Bucket3D.nY);
    const float stepZ= (Bucket3D.boxMax[2] - Bucket3D.boxMin[2]) / float(Bucket3D.nZ);
    // Set transformation
    glPushMatrix();
    glTranslatef(Bucket3D.boxMin[0] + 0.5f * stepX, Bucket3D.boxMin[1] + 0.5f * stepY, Bucket3D.boxMin[2] + 0.5f * stepZ);
    glScalef(stepX, stepY, stepZ);
    for (int x= 0; x < Bucket3D.nX; x++) {
      for (int y= 0; y < Bucket3D.nY; y++) {
        for (int z= 0; z < Bucket3D.nZ; z++) {
          // Color by occupancy
          float r= 0.5f, g= 0.5f, b= 0.5f;
          Colormap::RatioToJetBrightSmooth((float)Bucket3D.buckets[Bucket3D.GetFlatIdx(x, y, z)].size() / (float)Bucket3D.bucketCapacity, r, g, b);
          glColor3f(r, g, b);
          // Draw wire box
          glPushMatrix();
          glTranslatef((float)x, (float)y, (float)z);
          glutWireCube(0.99);
          glPopMatrix();
        }
      }
    }
    glPopMatrix();
    glLineWidth(1.0f);
  }


  // Draw triangles
  if (D.displayMode[4]) {
    if (!MetaballIsUpdated)
      ComputeMetaballs();
    glEnable(GL_LIGHTING);
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_TRIANGLES);
    for (int k= 0; k < (int)Tris.size(); k++) {
      Vec::Vec3 v0(Verts[Tris[k][0]][0], Verts[Tris[k][0]][1], Verts[Tris[k][0]][2]);
      Vec::Vec3 v1(Verts[Tris[k][1]][0], Verts[Tris[k][1]][1], Verts[Tris[k][1]][2]);
      Vec::Vec3 v2(Verts[Tris[k][2]][0], Verts[Tris[k][2]][1], Verts[Tris[k][2]][2]);
      Vec::Vec3 n0= (v1 - v0).cross(v2 - v0).normalized();
      Vec::Vec3 n1= n0;
      Vec::Vec3 n2= n0;
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
  else {
    MetaballIsUpdated= false;
    Verts.clear();
    Tris.clear();
  }

  // Write the status
  D.Status.clear();
  D.Status.resize(4);
  D.Status[0]= std::format("NbBuckets:{}x{}x{}={}", Bucket3D.nX, Bucket3D.nY, Bucket3D.nZ, Bucket3D.nXYZ);
  D.Status[1]= std::format("NbParticles:{}", (int)Pos.size());
  D.Status[2]= std::format("SimTime:{:.6f}ms", SimTime);
  if (BucketOverflown) D.Status[3]= std::string{"BUCKET OVERFLOW"};

  if (D.UI[VerboseLevel____].I() >= 1) printf("DrawT %f\n", Timer::PopTimer());
}
