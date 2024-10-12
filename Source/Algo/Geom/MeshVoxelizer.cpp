
#include "MeshVoxelizer.hpp"


// Standard lib
#include <array>
#include <cstdio>
#include <map>
#include <vector>

// Algo headers
#include "Geom/BoxGrid.hpp"
#include "Type/Field.hpp"
#include "Type/Vec.hpp"


void MeshVoxelizer::ComputeNarrowBandSignedDistanceField(const std::vector<std::array<double, 3>>& iVertices,
                                                         const std::vector<std::array<int, 3>>& iTriangles,
                                                         const std::array<double, 3>& iBoxMin,
                                                         const std::array<double, 3>& iBoxMax,
                                                         const int iNbX,
                                                         const int iNbY,
                                                         const int iNbZ,
                                                         const bool iIsCentered,
                                                         Field::Field3<char>& oIsSet,
                                                         Field::Field3<double>& oDist) {
  // Check inputs
  if (iVertices.size() < 3 || iTriangles.size() < 1 || iNbX <= 0 || iNbY <= 0 || iNbZ <= 0) {
    printf("[ERROR]: Invalid inputs for MeshVoxelizer::ComputeNarrowBandSignedDistanceField\n");
    return;
  }

  // Adjust the bounding box if computing voxel-centered distance values
  const double stepX= (iBoxMax[0] - iBoxMin[0]) / double(iNbX);
  const double stepY= (iBoxMax[1] - iBoxMin[1]) / double(iNbY);
  const double stepZ= (iBoxMax[2] - iBoxMin[2]) / double(iNbZ);
  const double bandDist= 1.01 * std::sqrt(stepX * stepX + stepY * stepY + stepZ * stepZ);
  std::array<double, 3> boxMin= iBoxMin;
  std::array<double, 3> boxMax= iBoxMax;
  if (iIsCentered) {
    boxMin= {iBoxMin[0] + 0.5 * stepX, iBoxMin[1] + 0.5 * stepY, iBoxMin[2] + 0.5 * stepZ};
    boxMax= {iBoxMax[0] - 0.5 * stepX, iBoxMax[1] - 0.5 * stepY, iBoxMax[2] - 0.5 * stepZ};
  }

  // Initialize variables
  oIsSet= Field::Field3<char>(iNbX, iNbY, iNbZ, 0);
  oDist= Field::Field3<double>(iNbX, iNbY, iNbZ, bandDist);
  Field::Field3<char> distSign(iNbX, iNbY, iNbZ, 0);

  // Compute face edge and vertex normals with incident angle weights
  // Needed to get correct sign in general 3D case
  // http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/1289/pdf/imm1289.pdf
  std::vector<Vec::Vec3<double>> faceNormals((int)iTriangles.size(), {0.0, 0.0, 0.0});
  std::vector<Vec::Vec3<double>> vertNormals((int)iVertices.size(), {0.0, 0.0, 0.0});
  std::map<std::pair<int, int>, Vec::Vec3<double>> edgeNormals;
  for (int k= 0; k < (int)iTriangles.size(); k++) {
    Vec::Vec3<double> AB= Vec::Vec3<double>(iVertices[iTriangles[k][1]]) - Vec::Vec3<double>(iVertices[iTriangles[k][0]]);
    Vec::Vec3<double> AC= Vec::Vec3<double>(iVertices[iTriangles[k][2]]) - Vec::Vec3<double>(iVertices[iTriangles[k][0]]);
    Vec::Vec3<double> BC= Vec::Vec3<double>(iVertices[iTriangles[k][2]]) - Vec::Vec3<double>(iVertices[iTriangles[k][1]]);
    faceNormals[k]= (AB).cross(AC).normalized();

    edgeNormals[{iTriangles[k][0], iTriangles[k][1]}].set(0.0, 0.0, 0.0);
    edgeNormals[{iTriangles[k][1], iTriangles[k][0]}].set(0.0, 0.0, 0.0);
    edgeNormals[{iTriangles[k][0], iTriangles[k][2]}].set(0.0, 0.0, 0.0);
    edgeNormals[{iTriangles[k][2], iTriangles[k][0]}].set(0.0, 0.0, 0.0);
    edgeNormals[{iTriangles[k][1], iTriangles[k][2]}].set(0.0, 0.0, 0.0);
    edgeNormals[{iTriangles[k][2], iTriangles[k][1]}].set(0.0, 0.0, 0.0);

    double angleA= std::acos(AB.dot(AC) / (AB.norm() * AC.norm()));
    double angleB= std::acos(BC.dot((-1.0 * AB)) / (BC.norm() * AB.norm()));
    double angleC= std::acos((-1.0 * AC).dot((-1.0 * BC)) / (AC.norm() * BC.norm()));
    vertNormals[iTriangles[k][0]]+= angleA * faceNormals[k];
    vertNormals[iTriangles[k][1]]+= angleB * faceNormals[k];
    vertNormals[iTriangles[k][2]]+= angleC * faceNormals[k];
  }
  for (int k= 0; k < (int)iTriangles.size(); k++) {
    edgeNormals[{iTriangles[k][0], iTriangles[k][1]}]+= faceNormals[k];
    edgeNormals[{iTriangles[k][1], iTriangles[k][0]}]+= faceNormals[k];
    edgeNormals[{iTriangles[k][0], iTriangles[k][2]}]+= faceNormals[k];
    edgeNormals[{iTriangles[k][2], iTriangles[k][0]}]+= faceNormals[k];
    edgeNormals[{iTriangles[k][1], iTriangles[k][2]}]+= faceNormals[k];
    edgeNormals[{iTriangles[k][2], iTriangles[k][1]}]+= faceNormals[k];
  }
  for (int k= 0; k < (int)iVertices.size(); k++)
    vertNormals[k].normalize();

  for (std::map<std::pair<int, int>, Vec::Vec3<double>>::iterator it= edgeNormals.begin(); it != edgeNormals.end(); it++)
    it->second.normalize();

  // Intersect the mesh with the voxel field
  for (int k= 0; k < (int)iTriangles.size(); k++) {
    // Find triangle bounding box augmented by iMaxDistance
    std::array<double, 3> posMin, posMax;
    posMin[0]= std::min(std::min(iVertices[iTriangles[k][0]][0], iVertices[iTriangles[k][1]][0]), iVertices[iTriangles[k][2]][0]) - bandDist;
    posMin[1]= std::min(std::min(iVertices[iTriangles[k][0]][1], iVertices[iTriangles[k][1]][1]), iVertices[iTriangles[k][2]][1]) - bandDist;
    posMin[2]= std::min(std::min(iVertices[iTriangles[k][0]][2], iVertices[iTriangles[k][1]][2]), iVertices[iTriangles[k][2]][2]) - bandDist;
    posMax[0]= std::max(std::max(iVertices[iTriangles[k][0]][0], iVertices[iTriangles[k][1]][0]), iVertices[iTriangles[k][2]][0]) + bandDist;
    posMax[1]= std::max(std::max(iVertices[iTriangles[k][0]][1], iVertices[iTriangles[k][1]][1]), iVertices[iTriangles[k][2]][1]) + bandDist;
    posMax[2]= std::max(std::max(iVertices[iTriangles[k][0]][2], iVertices[iTriangles[k][1]][2]), iVertices[iTriangles[k][2]][2]) + bandDist;

    // Find voxelized version of augmented bounding box
    std::array<int, 3> idxMin= {0, 0, 0}, idxMax= {0, 0, 0};
    if (boxMax[0] > boxMin[0]) idxMin[0]= std::max((int)std::round(double(iNbX - 1) * (posMin[0] - boxMin[0]) / (boxMax[0] - boxMin[0])), 0);
    if (boxMax[1] > boxMin[1]) idxMin[1]= std::max((int)std::round(double(iNbY - 1) * (posMin[1] - boxMin[1]) / (boxMax[1] - boxMin[1])), 0);
    if (boxMax[2] > boxMin[2]) idxMin[2]= std::max((int)std::round(double(iNbZ - 1) * (posMin[2] - boxMin[2]) / (boxMax[2] - boxMin[2])), 0);
    if (boxMax[0] > boxMin[0]) idxMax[0]= std::min((int)std::round(double(iNbX - 1) * (posMax[0] - boxMin[0]) / (boxMax[0] - boxMin[0])), iNbX - 1);
    if (boxMax[1] > boxMin[1]) idxMax[1]= std::min((int)std::round(double(iNbY - 1) * (posMax[1] - boxMin[1]) / (boxMax[1] - boxMin[1])), iNbY - 1);
    if (boxMax[2] > boxMin[2]) idxMax[2]= std::min((int)std::round(double(iNbZ - 1) * (posMax[2] - boxMin[2]) / (boxMax[2] - boxMin[2])), iNbZ - 1);

// Find neighboring voxels that intersect the current triangle
#pragma omp parallel for
    for (int x= idxMin[0]; x <= idxMax[0]; x++) {
      for (int y= idxMin[1]; y <= idxMax[1]; y++) {
        for (int z= idxMin[2]; z <= idxMax[2]; z++) {
          // Get current point position
          Vec::Vec3<double> cellPos(boxMin);
          if (iNbX > 1) cellPos[0]+= double(x) * (boxMax[0] - boxMin[0]) / double(iNbX - 1);
          if (iNbY > 1) cellPos[1]+= double(y) * (boxMax[1] - boxMin[1]) / double(iNbY - 1);
          if (iNbZ > 1) cellPos[2]+= double(z) * (boxMax[2] - boxMin[2]) / double(iNbZ - 1);

          // Get projection on triangle
          Vec::Vec3<double> projPoint;
          int elemType= -1;
          MeshVoxelizer::ComputePointToTriangleProjection(iVertices[iTriangles[k][0]], iVertices[iTriangles[k][1]], iVertices[iTriangles[k][2]], cellPos, projPoint, elemType);
          double unsignedDist= (cellPos - projPoint).norm();

          // Update distance field if needed
          if (unsignedDist <= bandDist) {
            Vec::Vec3<double> projNormal(0.0, 0.0, 0.0);
            if (elemType == 0) projNormal= faceNormals[k];
            if (elemType == 1) projNormal= vertNormals[iTriangles[k][0]];
            if (elemType == 2) projNormal= vertNormals[iTriangles[k][1]];
            if (elemType == 3) projNormal= vertNormals[iTriangles[k][2]];
            if (elemType == 4) projNormal= edgeNormals[{iTriangles[k][0], iTriangles[k][1]}];
            if (elemType == 5) projNormal= edgeNormals[{iTriangles[k][0], iTriangles[k][2]}];
            if (elemType == 6) projNormal= edgeNormals[{iTriangles[k][1], iTriangles[k][2]}];

            if (oDist.at(x, y, z) > unsignedDist) {
              oDist.at(x, y, z)= unsignedDist;
              oIsSet.at(x, y, z)= 1;
              distSign.at(x, y, z)= ((cellPos - projPoint).dot(projNormal) < 0.0) ? -1 : 1;
            }
          }
        }
      }
    }
  }

  // Apply the sign
  for (int xyz= 0; xyz < oDist.nXYZ; xyz++)
    if (distSign.at(xyz) < 0)
      oDist.at(xyz)= -oDist.at(xyz);
}


void MeshVoxelizer::ComputePointToTriangleProjection(const Vec::Vec3<double>& iTriangleA,
                                                     const Vec::Vec3<double>& iTriangleB,
                                                     const Vec::Vec3<double>& iTriangleC,
                                                     const Vec::Vec3<double>& iPoint,
                                                     Vec::Vec3<double>& oProjPoint,
                                                     int& oElemType) {
  // Implementation based on barycentric coordinates
  // http://web.mit.edu/ehliu/Public/Darmofal/DistancePoint3Triangle3.pdf
  // https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

  //      r2
  //
  //         C
  //             /|
  //            / | r1
  //   r3   t  /  |
  //          / r0|
  //         /____|
  //        A   s  B
  //
  //   r4      r5     r6
  //
  // REGION 0 : Proj on face   ABC: oElemType= 0
  // REGION 4 : Proj on vertex A:   oElemType= 1
  // REGION 6 : Proj on vertex B:   oElemType= 2
  // REGION 2 : Proj on vertex C:   oElemType= 3
  // REGION 5 : Proj on edge   AB:  oElemType= 4
  // REGION 3 : Proj on edge   AC:  oElemType= 5
  // REGION 1 : Proj on edge   BC:  oElemType= 6

  // Precomputations
  Vec::Vec3<double> edge0= iTriangleB - iTriangleA;
  Vec::Vec3<double> edge1= iTriangleC - iTriangleA;
  Vec::Vec3<double> v0= iPoint - iTriangleA;

  double a= edge0.dot(edge0);
  double b= edge0.dot(edge1);
  double c= edge1.dot(edge1);
  double d= -edge0.dot(v0);
  double e= -edge1.dot(v0);

  double s= b * e - c * d;
  double t= b * d - a * e;
  double det= a * c - b * b;

  // Barycentric coordinates to find closest point on triangle
  if (s + t <= det) {
    if (s < 0.0) {
      if (t < 0.0) {  // REGION 4 : Proj on vertex A:   oElemType= 1
        oElemType= 1;
        if (d < 0.0) {
          t= 0.0;
          if (-d >= a) s= 1.0;
          else s= -d / a;
        }
        else {
          s= 0.0;
          if (e >= 0.0) t= 0.0;
          else if (-e >= c) t= 1.0;
          else t= -e / c;
        }
      }
      else {  // REGION 3 : Proj on edge   AC:  oElemType= 5
        oElemType= 5;
        s= 0.0;
        if (e >= 0.0) t= 0.0;
        else if (-e >= c) t= 1.0;
        else t= -e / c;
      }
    }
    else if (t < 0.0) {  // REGION 5 : Proj on edge   AB:  oElemType= 4
      oElemType= 4;
      t= 0.0;
      if (d >= 0.0) s= 0.0;
      else if (-d >= a) s= 1.0;
      else s= -d / a;
    }
    else {  // REGION 0 : Proj on face   ABC: oElemType= 0
      oElemType= 0;
      s/= det;
      t/= det;
    }
  }
  else {
    if (s < 0.0) {  // REGION 2 : Proj on vertex C:   oElemType= 3
      oElemType= 3;
      double tmp0= b + d;
      double tmp1= c + e;
      if (tmp1 > tmp0) {
        double numer= tmp1 - tmp0;
        double denom= a - 2 * b + c;
        if (numer >= denom) {
          s= 1.0;
          t= 0.0;
        }
        else {
          s= numer / denom;
          t= 1.0 - s;
        }
      }
      else {
        s= 0.0;
        if (tmp1 <= 0.0) t= 1.0;
        else if (e >= 0) t= 0.0;
        else t= -e / c;
      }
    }
    else if (t < 0.0) {  // REGION 6 : Proj on vertex B:   oElemType= 2
      oElemType= 2;
      double tmp0= b + e;
      double tmp1= a + d;
      if (tmp1 > tmp0) {
        double numer= tmp1 - tmp0;
        double denom= a - 2.0 * b + c;
        if (numer >= denom) {
          t= 1.0;
          s= 0.0;
        }
        else {
          t= numer / denom;
          s= 1.0 - t;
        }
      }
      else {
        t= 0.0;
        if (tmp1 <= 0) s= 1.0;
        else if (d >= 0.0) s= 0.0;
        else s= -d / a;
      }
    }
    else {  // REGION 1 : Proj on edge   BC:  oElemType= 6
      oElemType= 6;
      double numer= c + e - b - d;
      if (numer <= 0.0) {
        s= 0.0;
        t= 1.0;
      }
      else {
        double denom= a - 2.0 * b + c;
        if (numer >= denom) {
          s= 1.0;
          t= 0.0;
        }
        else {
          s= numer / denom;
          t= 1.0 - s;
        }
      }
    }
  }

  // Compute coordinates of projected point
  oProjPoint= iTriangleA + s * edge0 + t * edge1;
}
