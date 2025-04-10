#include "ParticForceLaw.hpp"


// Standard lib
#include <cmath>
#include <vector>

// Algo headers
#include "FileIO/FileInput.hpp"
#include "Geom/Sketch.hpp"
#include "Type/Vec.hpp"
#include "Util/Timer.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void ParticForceLaw::BuildBaseCloud(std::vector<Vec::Vec3<float>>& oPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the base point cloud
  oPointCloud.clear();
  if (D.UI[LatticePitch____].F() <= 0.0f) return;
  std::vector<Vec::Vec3<float>> tmpPoints;

  // Regular cubic lattice patterns
  if (D.UI[LatticePattern__].I() == 0 || D.UI[LatticePattern__].I() == 1 || D.UI[LatticePattern__].I() == 2) {
    float minDist= 0.0f;
    if (D.UI[LatticePattern__].I() == 0) minDist= 1.0f;                                    // SCC pattern
    if (D.UI[LatticePattern__].I() == 1) minDist= (2.0f / 3.0f) * std::sqrt(3.0f) / 2.0f;  // BCC pattern
    if (D.UI[LatticePattern__].I() == 2) minDist= std::sqrt(2.0f) / 2.0f;                  // FCC pattern
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * minDist));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * minDist));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * minDist));
    for (int x= 0; x < (cellNbX / 2) * 2 + 1; x++) {
      for (int y= 0; y < (cellNbY / 2) * 2 + 1; y++) {
        for (int z= 0; z < (cellNbZ / 2) * 2 + 1; z++) {
          bool keep= false;
          if (D.UI[LatticePattern__].I() == 0) keep= true;                                                                                // SCC pattern
          if (D.UI[LatticePattern__].I() == 1 && ((x % 2 == 0 && y % 2 == 0 && z % 2 == 0) || (x % 2 + y % 2 + z % 2 == 3))) keep= true;  // BCC pattern
          if (D.UI[LatticePattern__].I() == 2 && ((x + y + z) % 2 == 0)) keep= true;                                                      // FCC pattern
          if (keep)
            tmpPoints.push_back(Vec::Vec3<float>(
                D.boxMin[0] + x * D.UI[LatticePitch____].F() * minDist,
                D.boxMin[1] + y * D.UI[LatticePitch____].F() * minDist,
                D.boxMin[2] + z * D.UI[LatticePitch____].F() * minDist));
        }
      }
    }
  }
  // HCP pattern with layers along X
  else if (D.UI[LatticePattern__].I() == 3) {
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * std::sqrt(6.0f) / 3.0f));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F()));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * 0.5f * std::sqrt(3.0f)));
    for (int x= 0; x < cellNbX + 1; x++)
      for (int y= 0; y < cellNbY; y++)
        for (int z= 0; z < cellNbZ; z++)
          tmpPoints.push_back(0.5f * D.UI[LatticePitch____].F() *
                                Vec::Vec3<float>(float(x) * 2.0f * std::sqrt(6.0f) / 3.0f,
                                                 2.0f * float(y) + float((z + x) % 2),
                                                 std::sqrt(3.0f) * (float(z) + float(x % 2) / 3.0f)));
  }
  // HCP pattern with layers along Y
  else if (D.UI[LatticePattern__].I() == 4) {
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F() * 0.5f * std::sqrt(3.0f)));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * std::sqrt(6.0f) / 3.0f));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F()));
    for (int x= 0; x < cellNbX; x++)
      for (int y= 0; y < cellNbY + 1; y++)
        for (int z= 0; z < cellNbZ; z++)
          tmpPoints.push_back(0.5f * D.UI[LatticePitch____].F() *
                                Vec::Vec3<float>(std::sqrt(3.0f) * (float(x) + float(y % 2) / 3.0f),
                                                 float(y) * 2.0f * std::sqrt(6.0f) / 3.0f,
                                                 2.0f * float(z) + float((x + y) % 2)));
  }
  // HCP pattern with layers along Z
  else if (D.UI[LatticePattern__].I() == 5) {
    const int cellNbX= (int)std::ceil(float(D.boxMax[0] - D.boxMin[0]) / (D.UI[LatticePitch____].F()));
    const int cellNbY= (int)std::ceil(float(D.boxMax[1] - D.boxMin[1]) / (D.UI[LatticePitch____].F() * 0.5f * std::sqrt(3.0f)));
    const int cellNbZ= (int)std::ceil(float(D.boxMax[2] - D.boxMin[2]) / (D.UI[LatticePitch____].F() * std::sqrt(6.0f) / 3.0f));
    for (int x= 0; x < cellNbX; x++)
      for (int y= 0; y < cellNbY; y++)
        for (int z= 0; z < cellNbZ + 1; z++)
          tmpPoints.push_back(0.5f * D.UI[LatticePitch____].F() *
                                Vec::Vec3<float>(2.0f * float(x) + float((y + z) % 2),
                                                 std::sqrt(3.0f) * (float(y) + float(z % 2) / 3.0f),
                                                 float(z) * 2.0f * std::sqrt(6.0f) / 3.0f));
  }

  // Recenter the point cloud then only keep the points in the bounding box
  Vec::Vec3<float> avgPos(0.0f, 0.0f, 0.0f);
  for (int k= 0; k < (int)tmpPoints.size(); k++)
    avgPos= avgPos + tmpPoints[k];
  avgPos/= (float)tmpPoints.size();
  for (int k= 0; k < (int)tmpPoints.size(); k++)
    for (int dim= 0; dim < 3; dim++)
      tmpPoints[k][dim]= tmpPoints[k][dim] - avgPos[dim] + (D.boxMin[dim] + D.boxMax[dim]) / 2.0f;
  for (int k= 0; k < (int)tmpPoints.size(); k++)
    if (tmpPoints[k][0] >= D.boxMin[0] && tmpPoints[k][1] >= D.boxMin[1] && tmpPoints[k][2] >= D.boxMin[2] &&
        tmpPoints[k][0] <= D.boxMax[0] && tmpPoints[k][1] <= D.boxMax[1] && tmpPoints[k][2] <= D.boxMax[2])
      oPointCloud.push_back(tmpPoints[k]);

  if (D.UI[VerboseLevel____].I() >= 1) printf("BaseCloudT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildScenario(const std::vector<Vec::Vec3<float>>& iPointCloud) {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Calculate some helper variables for scenario setup
  const Vec::Vec3<float> BoxMin(D.boxMin[0], D.boxMin[1], D.boxMin[2]);
  const Vec::Vec3<float> BoxMax(D.boxMax[0], D.boxMax[1], D.boxMax[2]);
  const float BoxDiag= (BoxMax - BoxMin).norm();
  const Vec::Vec3<float> BCForVecNega(-D.UI[BCForX__________].F(), -D.UI[BCForY__________].F(), -D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCForVecPosi(D.UI[BCForX__________].F(), D.UI[BCForY__________].F(), D.UI[BCForZ__________].F());
  const Vec::Vec3<float> BCVelVecNega(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F());
  const Vec::Vec3<float> BCVelVecPosi(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F());

  // Load the 2D scenario file if needed
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (D.UI[ScenarioPreset__].I() == 0) {
    if (D.UI[Scenario2DID____].I() == 0) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/Coupon.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 1) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/SandiaFracture.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 2) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/KalthoffFracture.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 3) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/Auxetic.bmp", imageRGBA, false);
    if (D.UI[Scenario2DID____].I() == 4) FileInput::LoadImageBMPFile("./FileInput/StrucScenarios/Logo.bmp", imageRGBA, false);
  }

  // Load the 3D scenario file if needed
  if (D.UI[ScenarioPreset__].I() == 1) {
    // TODO
  }

  // Add the subset of points for the current scenario
  for (int k= 0; k < (int)iPointCloud.size(); k++) {
    const Vec::Vec3<float> RelPos= (iPointCloud[k] - BoxMin).cwiseDiv(BoxMax - BoxMin);
    // 2D loaded scenario
    if (D.UI[ScenarioPreset__].I() == 0) {
      if (!imageRGBA.empty()) {
        if (std::abs(RelPos[0] - 0.5f) < 0.5f * D.UI[Scenario2DThick_].F()) {
          const float posW= (float)(imageRGBA.size() - 1) * RelPos[1];
          const float posH= (float)(imageRGBA[0].size() - 1) * RelPos[2];
          const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
          const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
          const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
          if (colRGBA[3] > 0.5f) {
            Pos.push_back(iPointCloud[k]);
            if (colRGBA[0] < 0.1f) Sensor.push_back(1);
            if (colRGBA[0] > 0.9f) BCPos.push_back(1);
            if (colRGBA[1] < 0.1f) BCVel.push_back(-1);
            if (colRGBA[1] > 0.9f) BCVel.push_back(1);
            if (colRGBA[2] < 0.1f) BCFor.push_back(-1);
            if (colRGBA[2] > 0.9f) BCFor.push_back(1);
          }
        }
      }
    }
    // 3D loaded scenario
    else if (D.UI[ScenarioPreset__].I() == 1) {
      // TODO set particle based on loaded mesh
    }
    // Full set of points
    else if (D.UI[ScenarioPreset__].I() == 2) {
      Pos.push_back(iPointCloud[k]);
    }
    // Box falling on steps
    else if (D.UI[ScenarioPreset__].I() == 3) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.8f)).abs().maxCoeff() < 0.1f * BoxDiag) {
        Pos.push_back(iPointCloud[k]);
        BCFor.push_back(-1);
      }
      else if (RelPos[0] > 0.5f && RelPos[2] > 0.3f && RelPos[2] < 0.5f && (RelPos[1] + RelPos[2]) <= 0.8f && RelPos[1] < 0.5f) {
        Pos.push_back(iPointCloud[k]);
        BCPos.push_back(1);
      }
      else if (RelPos[0] < 0.6f && RelPos[2] < 0.06f && RelPos[1] > 0.5f) {
        Pos.push_back(iPointCloud[k]);
        BCPos.push_back(1);
      }
    }
    // Ball blasting through wall
    else if (D.UI[ScenarioPreset__].I() == 4) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.2f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Mat.push_back(0);
      }
      else if (RelPos[0] > 0.1f && RelPos[0] < 0.9f &&
               RelPos[1] > 0.1f && RelPos[1] < 0.9f &&
               RelPos[2] > 0.45f && RelPos[2] < 0.55f) {
        Pos.push_back(iPointCloud[k]);
        Mat.push_back(1);
      }
    }
    // Balls colliding
    else if (D.UI[ScenarioPreset__].I() == 5) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.45f, 0.2f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Mat.push_back(0);
      }
      else if ((RelPos - Vec::Vec3<float>(0.5f, 0.55f, 0.8f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(-D.UI[BCVelX__________].F(), -D.UI[BCVelY__________].F(), -D.UI[BCVelZ__________].F()));
        Mat.push_back(1);
      }
    }
    // Coupon stretch - Velocity or Force driven
    else if (D.UI[ScenarioPreset__].I() == 6 || D.UI[ScenarioPreset__].I() == 7) {
      if (RelPos[0] > 0.1f && RelPos[0] < 0.9f) {
        if (RelPos[1] > 0.1f && RelPos[1] < 0.9f) {
          if (RelPos[2] > 0.2f && RelPos[2] < 0.8f) {
            Pos.push_back(iPointCloud[k]);
            if (D.UI[ScenarioPreset__].I() == 6) {
              if (RelPos[2] < 0.3f) {
                BCVel.push_back(-1);
              }
              else if (RelPos[2] > 0.7f) {
                Sensor.push_back(1);
                BCVel.push_back(1);
              }
            }
            else {
              if (RelPos[2] < 0.3f) {
                BCFor.push_back(-1);
              }
              else if (RelPos[2] > 0.7f) {
                Sensor.push_back(1);
                BCFor.push_back(1);
              }
            }
          }
        }
      }
    }
    // 3 Point flexture - Velocity or Force driven
    else if (D.UI[ScenarioPreset__].I() == 8 || D.UI[ScenarioPreset__].I() == 9) {
      if (RelPos[1] > 0.1f && RelPos[1] < 0.9f) {
        if ((RelPos - Vec::Vec3<float>(RelPos[0], 0.15f, 0.38f)).norm() < 0.05f ||
            (RelPos - Vec::Vec3<float>(RelPos[0], 0.85f, 0.38f)).norm() < 0.05f) {
          Pos.push_back(iPointCloud[k]);
          BCPos.push_back(1);
        }
        else if ((RelPos - Vec::Vec3<float>(RelPos[0], 0.5f, 0.62f)).norm() < 0.05f) {
          Pos.push_back(iPointCloud[k]);
          Sensor.push_back(1);
          if (D.UI[ScenarioPreset__].I() == 8) BCVel.push_back(-1);
          if (D.UI[ScenarioPreset__].I() == 9) BCFor.push_back(-1);
        }
        else if (RelPos[2] > 0.45f && RelPos[2] < 0.55f) {
          if (std::abs(RelPos[0] - 0.5f) < 0.5f * D.UI[Scenario2DThick_].F()) {
            Pos.push_back(iPointCloud[k]);
          }
        }
      }
    }
    // Simple sphere
    else if (D.UI[ScenarioPreset__].I() == 10) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.5f)).norm() < 0.4f) {
        Pos.push_back(iPointCloud[k]);
        BCFor.push_back(-1);
        Sensor.push_back(1);
      }
    }
    // Ball blasting through Bimaterial cylinder
    else if (D.UI[ScenarioPreset__].I() == 11) {
      if ((RelPos - Vec::Vec3<float>(0.5f, 0.5f, 0.15f)).norm() < 0.10f) {
        Pos.push_back(iPointCloud[k]);
        Vel.push_back(Vec::Vec3<float>(D.UI[BCVelX__________].F(), D.UI[BCVelY__________].F(), D.UI[BCVelZ__________].F()));
        Mat.push_back(0);
      }
      else if (RelPos[1] > 0.1f && RelPos[1] < 0.9f && (RelPos - Vec::Vec3<float>(0.5f, RelPos[1], 0.5f)).norm() < 0.05f) {
        Pos.push_back(iPointCloud[k]);
        Mat.push_back(0);
      }
      else if (RelPos[1] > 0.1f && RelPos[1] < 0.9f && (RelPos - Vec::Vec3<float>(0.5f, RelPos[1], 0.5f)).norm() < 0.15f) {
        Pos.push_back(iPointCloud[k]);
        Mat.push_back(1);
      }
    }

    // Fill the missing default values
    if (Mat.size() < Pos.size()) Mat.push_back(0);
    if (ForceMag.size() < Pos.size()) ForceMag.push_back(0.0);
    if (Neighbors.size() < Pos.size()) Neighbors.push_back(1);
    if (Sensor.size() < Pos.size()) Sensor.push_back(0);
    if (BCPos.size() < Pos.size()) BCPos.push_back(0);
    if (BCVel.size() < Pos.size()) BCVel.push_back(0);
    if (BCFor.size() < Pos.size()) BCFor.push_back(0);

    if (Ref.size() < Pos.size()) Ref.push_back(Pos[Pos.size() - 1]);
    if (Vel.size() < Pos.size()) Vel.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (For.size() < Pos.size()) For.push_back(Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
    if (Col.size() < Pos.size()) Col.push_back(Vec::Vec3<float>(0.5f, 0.5f, 0.5f));
  }

  SimTime= 0.0f;
  RunID++;

  if (D.UI[VerboseLevel____].I() >= 1) printf("ScenarioT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildForceLaws() {
  if (D.UI[VerboseLevel____].I() >= 1) Timer::PushTimer();

  // Reset the force laws
  constexpr int nbForceLaws= 2;
  ForceLaw.clear();
  ForceLaw.resize(nbForceLaws);
  ForceLawStep.clear();
  ForceLawStep.resize(nbForceLaws);
  ForceLawRange.clear();
  ForceLawRange.resize(nbForceLaws);

  // Create the force laws
  for (int idxMat= 0; idxMat < (int)ForceLaw.size(); idxMat++) {
    int presetID= -1;
    float forceLawScale= 1.0f;
    if (idxMat == 0) {
      presetID= D.UI[ForceLawPresetA_].I();
      forceLawScale= D.UI[ForceLawScaleA__].F();
    }
    if (idxMat == 1) {
      presetID= D.UI[ForceLawPresetB_].I();
      forceLawScale= D.UI[ForceLawScaleB__].F();
    }
    // Custom force law
    if (presetID == 0) {
      BuildForceLawPolyline(D.UI[ForceLawA_0_00__].D(), D.UI[ForceLawA_0_80__].D(), D.UI[ForceLawA_0_90__].D(), D.UI[ForceLawA_0_95__].D(),
                            D.UI[ForceLawA_1_00__].D(), D.UI[ForceLawA_1_05__].D(), D.UI[ForceLawA_1_10__].D(), D.UI[ForceLawA_1_20__].D(),
                            D.UI[ForceLawA_1_30__].D(), D.UI[ForceLawA_1_40__].D(), D.UI[ForceLawA_1_50__].D(), D.UI[ForceLawA_2_00__].D(),
                            D.UI[ForceLawA_2_50__].D(), D.UI[ForceLawA_3_00__].D(), idxMat);
    }
    else if (presetID == 1) {
      BuildForceLawPolyline(D.UI[ForceLawB_0_00__].D(), D.UI[ForceLawB_0_80__].D(), D.UI[ForceLawB_0_90__].D(), D.UI[ForceLawB_0_95__].D(),
                            D.UI[ForceLawB_1_00__].D(), D.UI[ForceLawB_1_05__].D(), D.UI[ForceLawB_1_10__].D(), D.UI[ForceLawB_1_20__].D(),
                            D.UI[ForceLawB_1_30__].D(), D.UI[ForceLawB_1_40__].D(), D.UI[ForceLawB_1_50__].D(), D.UI[ForceLawB_2_00__].D(),
                            D.UI[ForceLawB_2_50__].D(), D.UI[ForceLawB_3_00__].D(), idxMat);
    }
    // Hard coded force laws from MIT paper
    else {
      ForceLawStep[idxMat]= 0.05f;
      ForceLaw[idxMat]= std::vector<float>{1.0f};
      //                                                      0.00,      0.05,      0.10,      0.15,      0.20,      0.25,      0.30,      0.35,      0.40,      0.45,      0.50,      0.55,      0.60,      0.65,      0.70,      0.75,      0.80,      0.85,      0.90,      0.95,      1.00,      1.05,      1.10,      1.15,      1.20,      1.25,      1.30,      1.35,      1.40,      1.45,      1.50,      1.55,      1.60,      1.65,      1.70,      1.75,      1.80,      1.85,      1.90,      1.95,      2.00,      2.05,      2.10,      2.15,      2.20,      2.25,      2.30,      2.35,      2.40,      2.45,      2.50,      2.55,      2.60,      2.65,      2.70,      2.75,      2.80,      2.85,      2.90,      2.95,      3.00
      if (presetID == 2) ForceLaw[idxMat]= std::vector<float>{1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +5.00E+05, +0.00E+00, -5.00E+05, -1.00E+06, -5.00E+05, +0.00E+00, +5.00E+05, +1.00E+06, +5.00E+05, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Elastic material
      if (presetID == 3) ForceLaw[idxMat]= std::vector<float>{1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +1.00E+06, +5.00E+05, +0.00E+00, -1.00E+05, +0.00E+00, +7.00E+04, +1.00E+05, +1.20E+05, +8.00E+04, +3.00E+04, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Brittle Material
      if (presetID == 4) ForceLaw[idxMat]= std::vector<float>{1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +1.00E+02, +9.50E+01, +8.80E+01, +7.00E+01, +5.00E+01, +2.50E+01, +0.00E+00, -2.00E+01, -2.80E+01, -3.00E+01, -2.80E+01, -2.50E+01, -2.00E+01, -1.20E+01, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Viscous material
      if (presetID == 5) ForceLaw[idxMat]= std::vector<float>{3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +3.00E+10, +2.00E+10, +0.00E+00, -2.00E+10, -2.50E+10, -2.00E+10, +0.00E+00, +2.00E+10, +2.50E+10, +2.00E+10, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Steel AISI 4340
      if (presetID == 6) ForceLaw[idxMat]= std::vector<float>{1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.50E+08, +1.00E+08, +0.00E+00, -1.50E+08, -2.30E+08, -1.50E+08, +0.00E+00, +1.10E+08, +1.00E+08, +5.00E+07, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00, +0.00E+00};  // Sample force law
      if (presetID == 7) ForceLaw[idxMat]= std::vector<float>{8.80E+04, +8.80E+04, +8.80E+04, +8.70E+04, +8.70E+04, +8.60E+04, +8.50E+04, +8.50E+04, +8.30E+04, +8.10E+04, +7.80E+04, +7.50E+04, +7.20E+04, +6.60E+04, +6.00E+04, +5.30E+04, +4.50E+04, +3.70E+04, +2.90E+04, +2.00E+04, +0.00E+00, -3.00E+04, -3.20E+04, -3.30E+04, -3.30E+04, -3.30E+04, -3.20E+04, -3.20E+04, -3.10E+04, -3.10E+04, -3.00E+04, -3.00E+04, -3.00E+04, -2.90E+04, -2.90E+04, -2.80E+04, -2.80E+04, -2.80E+04, -2.80E+04, -2.90E+04, -2.90E+04, -3.00E+04, -3.00E+04, -3.00E+04, -3.10E+04, -3.20E+04, -3.30E+04, -3.50E+04, -3.60E+04, -3.80E+04, -4.00E+04, -4.20E+04, -4.30E+04, -4.40E+04, -4.40E+04, -4.30E+04, -4.10E+04, -3.80E+04, -3.20E+04, -2.10E+04, +0.00E+00};  // Delrin force law optimized on force displacement curve
    }

    // Compute the effective range of the force law ignoring the zero tail
    int tailStartIdx= 0;
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      if (std::abs(ForceLaw[idxMat][k]) > 0.0f)
        tailStartIdx= k;
    ForceLawRange[idxMat]= ForceLawStep[idxMat] * (float)(tailStartIdx + 1) * D.UI[LatticePitch____].F();

    // Normalize the force law by the first value
    const float BaseVal= ForceLaw[idxMat][0];
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      ForceLaw[idxMat][k]/= BaseVal;

    // Scale the force law by the UI coeff
    for (int k= 0; k < (int)ForceLaw[idxMat].size(); k++)
      ForceLaw[idxMat][k]*= forceLawScale;
  }

  if (D.UI[VerboseLevel____].I() >= 1) printf("ForceLawsT %f\n", Timer::PopTimer());
}


void ParticForceLaw::BuildForceLawPolyline(const double v0_00, const double v0_80, const double v0_90, const double v0_95,
                                           const double v1_00, const double v1_05, const double v1_10, const double v1_20,
                                           const double v1_30, const double v1_40, const double v1_50, const double v2_00,
                                           const double v2_50, const double v3_00, const int iIdxMat) {
  // Create the force law as a smoothed polyline
  std::vector<std::array<double, 3>> PolylineA;
  PolylineA.push_back(std::array<double, 3>{0.0, 0.00, v0_00});
  PolylineA.push_back(std::array<double, 3>{0.0, 0.80, v0_80});
  PolylineA.push_back(std::array<double, 3>{0.0, 0.90, v0_90});
  PolylineA.push_back(std::array<double, 3>{0.0, 0.95, v0_95});
  PolylineA.push_back(std::array<double, 3>{0.0, 1.00, v1_00});
  std::vector<std::array<double, 3>> PolylineB;
  PolylineB.push_back(std::array<double, 3>{0.0, 1.05, v1_05});
  PolylineB.push_back(std::array<double, 3>{0.0, 1.10, v1_10});
  PolylineB.push_back(std::array<double, 3>{0.0, 1.20, v1_20});
  PolylineB.push_back(std::array<double, 3>{0.0, 1.30, v1_30});
  PolylineB.push_back(std::array<double, 3>{0.0, std::sqrt(2), v1_40});
  if (v1_40 != 0.0 || v1_50 != 0.0 || v2_00 != 0.0 || v2_50 != 0.0 || v3_00 != 0.0) {
    PolylineB.push_back(std::array<double, 3>{0.0, 1.50, v1_50});
    PolylineB.push_back(std::array<double, 3>{0.0, 2.00, v2_00});
    PolylineB.push_back(std::array<double, 3>{0.0, 2.50, v2_50});
    PolylineB.push_back(std::array<double, 3>{0.0, 3.00, v3_00});
  }
  Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineA);
  Sketch::PolylineSubdivideAndSmooth(true, 5, 5, PolylineB);
  std::vector<std::array<double, 3>> Polyline;
  Polyline.insert(Polyline.end(), PolylineA.begin(), PolylineA.end());
  Polyline.insert(Polyline.end(), PolylineB.begin(), PolylineB.end());

  // Sample the force value at fixed distance intervals
  constexpr int nbSamples= 211;  // 3*70+1 chosen to have a sample point at 1.0 and another very close to sqrt(2)
  constexpr float maxReach= 3.0f;
  ForceLaw[iIdxMat].resize(nbSamples);
  ForceLawStep[iIdxMat]= maxReach / float(nbSamples - 1);
  for (int k= 0; k < nbSamples; k++) {
    float dist= ForceLawStep[iIdxMat] * float(k);
    int idxLow= 0, idxUpp= 0;
    for (int idxVert= 0; idxVert < (int)Polyline.size() - 1; idxVert++) {
      idxLow= idxVert;
      idxUpp= std::min(idxVert + 1, (int)Polyline.size() - 1);
      if (Polyline[idxLow][1] <= dist && Polyline[idxUpp][1] >= dist) break;
    }
    double ratio= 0.0;
    if (Polyline[idxUpp][1] - Polyline[idxLow][1] > 0.0)
      ratio= (dist - Polyline[idxLow][1]) / (Polyline[idxUpp][1] - Polyline[idxLow][1]);
    ratio= std::min(std::max(ratio, 0.0), 1.0);
    ForceLaw[iIdxMat][k]= (1.0 - ratio) * Polyline[idxLow][2] + ratio * Polyline[idxUpp][2];
  }
}
