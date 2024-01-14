
#include "PrimitiveCSG.hpp"

// Standard lib
#include <array>
#include <cmath>
#include <vector>

// Eigen lib
#include <Eigen/Core>
#include <Eigen/Dense>

// Algo headers
#include "Geom/BoxGrid.hpp"
#include "Math/Field.hpp"


void PrimitiveCSG::Sphere(
    std::array<double, 3> const& iCenter,
    double const& iRadius,
    PrimitiveCSG::BooleanMode const& iMode,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  double stepX, stepY, stepZ, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  // Get primitives properties
  Eigen::Vector3d center(iCenter[0], iCenter[1], iCenter[2]);
  // Sweep through the field
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        // Get current point coordinates
        Eigen::Vector3d P(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
        // Compute the distance
        double distVal= (P - center).norm() - iRadius;
        // Update the distance
        if (iMode == PrimitiveCSG::BooleanMode::Union)
          ioDistanceField[x][y][z]= std::min(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Difference)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], -distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
          ioDistanceField[x][y][z]= std::max(-ioDistanceField[x][y][z], distVal);
      }
    }
  }
}


void PrimitiveCSG::Cylinder(
    std::array<double, 3> const& iPointA,
    std::array<double, 3> const& iPointB,
    double const& iRadiusA,
    double const& iRadiusB,
    bool const& iUseRoundedEnds,
    PrimitiveCSG::BooleanMode const& iMode,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  double stepX, stepY, stepZ, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  // Get primitives properties
  Eigen::Vector3d pointA(iPointA[0], iPointA[1], iPointA[2]);
  Eigen::Vector3d pointB(iPointB[0], iPointB[1], iPointB[2]);
  // Sweep through the field
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        // Get current point coordinates
        Eigen::Vector3d P(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
        // Compute the distance
        double distVal;
        if (iUseRoundedEnds) {
          distVal= std::min((pointA - P).norm() - iRadiusA, (pointB - P).norm() - iRadiusB);
          double relativePos= (P - pointA).dot(pointB - pointA) / (pointB - pointA).squaredNorm();
          if (relativePos > 0.0 && relativePos < 1.0) {
            Eigen::Vector3d projOnAxis= pointA + relativePos * (pointB - pointA);
            distVal= std::min(distVal, (P - projOnAxis).norm() - ((1.0 - relativePos) * iRadiusA + relativePos * iRadiusB));
          }
        }
        else {
          double relativePos= (P - pointA).dot(pointB - pointA) / (pointB - pointA).squaredNorm();
          Eigen::Vector3d projOnAxis= pointA + relativePos * (pointB - pointA);
          distVal= (P - projOnAxis).norm() - ((1.0 - relativePos) * iRadiusA + relativePos * iRadiusB);
          distVal= std::max(distVal, (relativePos < 0.0) ? (pointA - projOnAxis).norm() : -(pointA - projOnAxis).norm());
          distVal= std::max(distVal, (relativePos > 1.0) ? (pointB - projOnAxis).norm() : -(pointB - projOnAxis).norm());
        }
        // Update the distance
        if (iMode == PrimitiveCSG::BooleanMode::Union)
          ioDistanceField[x][y][z]= std::min(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Difference)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], -distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
          ioDistanceField[x][y][z]= std::max(-ioDistanceField[x][y][z], distVal);
      }
    }
  }
}


void PrimitiveCSG::ConeRound(
    std::array<double, 3> const& iPointA,
    std::array<double, 3> const& iPointB,
    double const& iRadiusA,
    double const& iRadiusB,
    PrimitiveCSG::BooleanMode const& iMode,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  double stepX, stepY, stepZ, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  // Get primitives properties
  Eigen::Vector3d pointA(iPointA[0], iPointA[1], iPointA[2]);
  Eigen::Vector3d pointB(iPointB[0], iPointB[1], iPointB[2]);
  // Sampling independent computations
  // Adapted from https://iquilezles.org/articles/distfunctions/
  Eigen::Vector3d ba= pointB - pointA;
  double l2= ba.dot(ba);
  double rr= iRadiusA - iRadiusB;
  double a2= l2 - rr * rr;
  double il2= 1.0 / l2;
  // Sweep through the field
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        // Get current point coordinates
        Eigen::Vector3d P(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
        // Sampling dependent computations
        Eigen::Vector3d pa= P - pointA;
        double y1= pa.dot(ba);
        double z1= y1 - l2;
        double x2= (pa * l2 - ba * y1).dot(pa * l2 - ba * y1);
        double y2= y1 * y1 * l2;
        double z2= z1 * z1 * l2;
        // Compute the distance
        double distVal;
        double k= copysign(1.0, rr) * rr * rr * x2;
        if (copysign(1.0, z1) * a2 * z2 > k)
          distVal= sqrt(x2 + z2) * il2 - iRadiusB;
        else if (copysign(1.0, y1) * a2 * y2 < k)
          distVal= sqrt(x2 + y2) * il2 - iRadiusA;
        else
          distVal= (sqrt(x2 * a2 * il2) + y1 * rr) * il2 - iRadiusA;
        // Update the distance
        if (iMode == PrimitiveCSG::BooleanMode::Union)
          ioDistanceField[x][y][z]= std::min(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Difference)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], -distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
          ioDistanceField[x][y][z]= std::max(-ioDistanceField[x][y][z], distVal);
      }
    }
  }
}


void PrimitiveCSG::AxisAlignedBox(
    std::array<double, 3> const& iCornerMin,
    std::array<double, 3> const& iCornerMax,
    PrimitiveCSG::BooleanMode const& iMode,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  double stepX, stepY, stepZ, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  // Get primitives properties
  Eigen::Vector3d cornerMin(iCornerMin[0], iCornerMin[1], iCornerMin[2]);
  Eigen::Vector3d cornerMax(iCornerMax[0], iCornerMax[1], iCornerMax[2]);
  // Sweep through the field
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        // Get current point coordinates
        Eigen::Vector3d P(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
        // Compute the distance
        Eigen::Vector3d d= (cornerMin - P).cwiseMax(P - cornerMax);
        double distVal= (Eigen::Vector3d(0.0, 0.0, 0.0).cwiseMax(d)).norm() + std::min(0.0, d.maxCoeff());
        // Update the distance
        if (iMode == PrimitiveCSG::BooleanMode::Union)
          ioDistanceField[x][y][z]= std::min(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Difference)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], -distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
          ioDistanceField[x][y][z]= std::max(-ioDistanceField[x][y][z], distVal);
      }
    }
  }
}


void PrimitiveCSG::AxisAlignedExtrudedSketch(
    int const iExtrusionDir,
    double const iExtruLimitMin,
    double const iExtruLimitMax,
    std::vector<std::array<double, 3>> const& iVertices,
    PrimitiveCSG::BooleanMode const& iMode,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Check inputs
  if (int(iVertices.size()) < 3) return;
  if (iExtrusionDir != 0 && iExtrusionDir != 1 && iExtrusionDir != 2) return;
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  double stepX, stepY, stepZ, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  // Build the 2D polygon in the sketch plane
  std::vector<Eigen::Vector2d> vertices2D(iVertices.size());
  for (int k= 0; k < (int)iVertices.size(); k++) {
    if (iExtrusionDir == 0) vertices2D[k]= {iVertices[k][1], iVertices[k][2]};
    else if (iExtrusionDir == 1) vertices2D[k]= {iVertices[k][0], iVertices[k][2]};
    else if (iExtrusionDir == 2) vertices2D[k]= {iVertices[k][0], iVertices[k][1]};
  }
  // Sweep through the sketch plane
  for (int x= 0; x < ((iExtrusionDir == 0) ? 1 : nbX); x++) {
    for (int y= 0; y < ((iExtrusionDir == 1) ? 1 : nbY); y++) {
      for (int z= 0; z < ((iExtrusionDir == 2) ? 1 : nbZ); z++) {
        // Get position in the sketch plane
        Eigen::Vector2d P;
        if (iExtrusionDir == 0) P= Eigen::Vector2d(double(y) * stepY + startY, double(z) * stepZ + startZ);
        else if (iExtrusionDir == 1) P= Eigen::Vector2d(double(x) * stepX + startX, double(z) * stepZ + startZ);
        else if (iExtrusionDir == 2) P= Eigen::Vector2d(double(x) * stepX + startX, double(y) * stepY + startY);
        // Compute the distance in the sketch plane
        // Adapted from https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
        double d= (P - vertices2D[0]).dot(P - vertices2D[0]);
        bool isPositive= true;
        for (int i= 0, j= int(vertices2D.size()) - 1; i < int(vertices2D.size()); j= i, i++) {
          Eigen::Vector2d e= vertices2D[j] - vertices2D[i];
          Eigen::Vector2d w= P - vertices2D[i];
          Eigen::Vector2d b= w - e * std::min(std::max(w.dot(e) / e.dot(e), 0.0), 1.0);
          d= std::min(d, b.dot(b));
          const bool condA= P[1] >= vertices2D[i][1];
          const bool condB= P[1] < vertices2D[j][1];
          const bool condC= e[0] * w[1] > e[1] * w[0];
          if ((condA && condB && condC) || (!condA && !condB && !condC))
            isPositive= !isPositive;
        }
        double sketchDistVal= std::sqrt(d);
        if (!isPositive)
          sketchDistVal= -sketchDistVal;
        // Sweep through the extrusion direction
        for (int xExt= ((iExtrusionDir == 0) ? 0 : x); xExt < ((iExtrusionDir == 0) ? nbX : x + 1); xExt++) {
          for (int yExt= ((iExtrusionDir == 1) ? 0 : y); yExt < ((iExtrusionDir == 1) ? nbY : y + 1); yExt++) {
            for (int zExt= ((iExtrusionDir == 2) ? 0 : z); zExt < ((iExtrusionDir == 2) ? nbZ : z + 1); zExt++) {
              // Limit the extrusion between the min and max limits
              double extruPos;
              if (iExtrusionDir == 0) extruPos= double(xExt) * stepX + startX;
              else if (iExtrusionDir == 1) extruPos= double(yExt) * stepY + startY;
              else if (iExtrusionDir == 2) extruPos= double(zExt) * stepZ + startZ;
              double extruDistVal= std::max(iExtruLimitMin - extruPos, extruPos - iExtruLimitMax);
              // Compute the distance
              double distVal= std::max(sketchDistVal, extruDistVal);
              // Update the distance
              if (iMode == PrimitiveCSG::BooleanMode::Union)
                ioDistanceField[xExt][yExt][zExt]= std::min(ioDistanceField[xExt][yExt][zExt], distVal);
              else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
                ioDistanceField[xExt][yExt][zExt]= std::max(ioDistanceField[xExt][yExt][zExt], distVal);
              else if (iMode == PrimitiveCSG::BooleanMode::Difference)
                ioDistanceField[xExt][yExt][zExt]= std::max(ioDistanceField[xExt][yExt][zExt], -distVal);
              else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
                ioDistanceField[xExt][yExt][zExt]= std::max(-ioDistanceField[xExt][yExt][zExt], distVal);
            }
          }
        }
      }
    }
  }
}


void PrimitiveCSG::GeneralExtrudedSketch(
    std::array<double, 3> const& iSketchCenter,
    std::array<double, 3> const& iSketchVector,
    double const iExtruLimitNega,
    double const iExtruLimitPosi,
    std::vector<std::array<double, 3>> const& iVertices,
    PrimitiveCSG::BooleanMode const& iMode,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Check inputs
  if (int(iVertices.size()) < 3) return;
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  double stepX, stepY, stepZ, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, true, startX, startY, startZ);
  Eigen::Vector3d sketchCenter(iSketchCenter[0], iSketchCenter[1], iSketchCenter[2]);
  Eigen::Vector3d sketchVector(iSketchVector[0], iSketchVector[1], iSketchVector[2]);
  // Compute rotation matrix
  // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
  Eigen::Matrix3d RotMat;
  {
    Eigen::Vector3d vertical(0.0, 0.0, 1.0);
    Eigen::Vector3d v= sketchVector.normalized().cross(vertical);
    double s= v.norm();
    double c= sketchVector.normalized().dot(vertical);
    Eigen::Matrix3d vskew;
    vskew << 0.0, -v[2], v[1],
        v[2], 0.0, -v[0],
        -v[1], v[0], 0.0;
    RotMat= Eigen::Matrix3d::Identity() + vskew + ((1 - c) / (s * s)) * vskew * vskew;
  }
  std::vector<Eigen::Vector2d> sketchVertices;
  for (int k= 0; k < int(iVertices.size()); k++) {
    Eigen::Vector3d pos3D(iVertices[k][0], iVertices[k][1], iVertices[k][2]);
    Eigen::Vector3d pos3DTranslat= pos3D - sketchCenter;
    Eigen::Vector3d pos3DRot= RotMat * pos3DTranslat;
    Eigen::Vector2d pos2D(pos3DRot[0], pos3DRot[1]);
    sketchVertices.push_back(pos2D);
  }
  // Sweep through the field
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        // Get position in the sketch plane
        Eigen::Vector3d pos3D(double(x) * stepX + startX, double(y) * stepY + startY, double(z) * stepZ + startZ);
        Eigen::Vector3d pos3DTranslat= pos3D - sketchCenter;
        Eigen::Vector3d pos3DRot= RotMat * pos3DTranslat;
        Eigen::Vector2d pos2D(pos3DRot[0], pos3DRot[1]);
        // Compute the distance in the sktech plane
        // https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
        double d= (pos2D - sketchVertices[0]).dot(pos2D - sketchVertices[0]);
        bool isPositive= true;
        for (int i= 0, j= int(sketchVertices.size()) - 1; i < int(sketchVertices.size()); j= i, i++) {
          Eigen::Vector2d e= sketchVertices[j] - sketchVertices[i];
          Eigen::Vector2d w= pos2D - sketchVertices[i];
          Eigen::Vector2d b= w - e * std::min(std::max(w.dot(e) / e.dot(e), 0.0), 1.0);
          d= std::min(d, b.dot(b));
          bool condA= pos2D[1] >= sketchVertices[i][1];
          bool condB= pos2D[1] < sketchVertices[j][1];
          bool condC= e[0] * w[1] > e[1] * w[0];
          if ((condA && condB && condC) || (!condA && !condB && !condC))
            isPositive= !isPositive;
        }
        double sketchDistVal= std::sqrt(d);
        if (!isPositive)
          sketchDistVal= -sketchDistVal;
        // Compute the distance orthogonal to the sketch plane
        double extruDistVal= std::max(iExtruLimitNega - pos3DRot[2], pos3DRot[2] - iExtruLimitPosi);
        // Limit the extrusion between the min and max limits
        double distVal= std::max(sketchDistVal, extruDistVal);
        // Update the distance
        if (iMode == PrimitiveCSG::BooleanMode::Union)
          ioDistanceField[x][y][z]= std::min(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::Difference)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], -distVal);
        else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
          ioDistanceField[x][y][z]= std::max(-ioDistanceField[x][y][z], distVal);
      }
    }
  }
}


void PrimitiveCSG::BooleanOperation(
    PrimitiveCSG::BooleanMode const& iMode,
    std::vector<std::vector<std::vector<double>>> const& iDistanceField,
    std::vector<std::vector<std::vector<double>>>& ioDistanceField) {
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetFieldDimensions(ioDistanceField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return;
  // Sweep through the field
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        // Update the distance
        if (iMode == PrimitiveCSG::BooleanMode::Union)
          ioDistanceField[x][y][z]= std::min(ioDistanceField[x][y][z], iDistanceField[x][y][z]);
        else if (iMode == PrimitiveCSG::BooleanMode::Intersection)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], iDistanceField[x][y][z]);
        else if (iMode == PrimitiveCSG::BooleanMode::Difference)
          ioDistanceField[x][y][z]= std::max(ioDistanceField[x][y][z], -iDistanceField[x][y][z]);
        else if (iMode == PrimitiveCSG::BooleanMode::DifferenceFlipOrder)
          ioDistanceField[x][y][z]= std::max(-ioDistanceField[x][y][z], iDistanceField[x][y][z]);
      }
    }
  }
}
