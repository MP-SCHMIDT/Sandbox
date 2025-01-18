
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

  // Precompute dimensions and box values
  double stepX, stepY, stepZ, stepDiag, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(iNbX, iNbY, iNbZ, iBoxMin, iBoxMax, iIsCentered, stepX, stepY, stepZ, stepDiag);
  BoxGrid::GetVoxelStart(iBoxMin, stepX, stepY, stepZ, iIsCentered, startX, startY, startZ);
  const double bandDist= 1.01 * stepDiag;
  std::array<int, 3> nb({iNbX, iNbY, iNbZ});
  std::array<double, 3> start({startX, startY, startZ});
  std::array<double, 3> gridBoxSizeInv({1.0/(double(iNbX-1)*stepX), 1.0/(double(iNbY-1)*stepY), 1.0/(double(iNbZ-1)*stepZ)});

  // Initialize variables
  oIsSet= Field::Field3<char>(iNbX, iNbY, iNbZ, 0);
  oDist= Field::Field3<double>(iNbX, iNbY, iNbZ, std::numeric_limits<double>::max());
  Field::Field3<char> distSign(iNbX, iNbY, iNbZ, 0);

  // Compute face edge and vertex normals with incident angle weights
  // Needed to get correct sign in general 3D case
  // http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/1289/pdf/imm1289.pdf
  std::vector<Vec::Vec3<double>> faceNormals((int)iTriangles.size(), {0.0, 0.0, 0.0});
  std::vector<Vec::Vec3<double>> vertNormals((int)iVertices.size(), {0.0, 0.0, 0.0});
  std::map<std::pair<int, int>, Vec::Vec3<double>> edgeNormals;
  for (int k= 0; k < (int)iTriangles.size(); k++) {
    const std::array<int, 3> tri= iTriangles[k];
    const Vec::Vec3<double> AB= Vec::Vec3<double>(iVertices[tri[1]]) - Vec::Vec3<double>(iVertices[tri[0]]);
    const Vec::Vec3<double> AC= Vec::Vec3<double>(iVertices[tri[2]]) - Vec::Vec3<double>(iVertices[tri[0]]);
    const Vec::Vec3<double> BC= Vec::Vec3<double>(iVertices[tri[2]]) - Vec::Vec3<double>(iVertices[tri[1]]);
    faceNormals[k]= (AB).cross(AC).normalized();

    edgeNormals[{tri[0], tri[1]}].set(0.0, 0.0, 0.0);
    edgeNormals[{tri[1], tri[0]}].set(0.0, 0.0, 0.0);
    edgeNormals[{tri[0], tri[2]}].set(0.0, 0.0, 0.0);
    edgeNormals[{tri[2], tri[0]}].set(0.0, 0.0, 0.0);
    edgeNormals[{tri[1], tri[2]}].set(0.0, 0.0, 0.0);
    edgeNormals[{tri[2], tri[1]}].set(0.0, 0.0, 0.0);

    const double angleA= std::acos(AB.dot(AC) / (AB.norm() * AC.norm()));
    const double angleB= std::acos(BC.dot((-1.0 * AB)) / (BC.norm() * AB.norm()));
    const double angleC= std::acos((-1.0 * AC).dot((-1.0 * BC)) / (AC.norm() * BC.norm()));
    vertNormals[tri[0]]+= angleA * faceNormals[k];
    vertNormals[tri[1]]+= angleB * faceNormals[k];
    vertNormals[tri[2]]+= angleC * faceNormals[k];
  }
  for (int k= 0; k < (int)iTriangles.size(); k++) {
    const std::array<int, 3> tri= iTriangles[k];
    edgeNormals[{tri[0], tri[1]}]+= faceNormals[k];
    edgeNormals[{tri[1], tri[0]}]+= faceNormals[k];
    edgeNormals[{tri[0], tri[2]}]+= faceNormals[k];
    edgeNormals[{tri[2], tri[0]}]+= faceNormals[k];
    edgeNormals[{tri[1], tri[2]}]+= faceNormals[k];
    edgeNormals[{tri[2], tri[1]}]+= faceNormals[k];
  }
  for (int k= 0; k < (int)iVertices.size(); k++)
    vertNormals[k].normalize();

  for (std::map<std::pair<int, int>, Vec::Vec3<double>>::iterator it= edgeNormals.begin(); it != edgeNormals.end(); it++)
    it->second.normalize();

  // Intersect the mesh with the voxel field
  for (unsigned int k= 0; k < iTriangles.size(); k++) {
    const std::array<int, 3> tri= iTriangles[k];
    // Find triangle bounding box augmented by the narrow band distance
    std::array<double, 3> posMin, posMax;
    for (unsigned int dim= 0; dim < 3; dim++) {
      posMin[dim]= std::min(std::min(iVertices[tri[0]][dim], iVertices[tri[1]][dim]), iVertices[tri[2]][dim]) - bandDist;
      posMax[dim]= std::max(std::max(iVertices[tri[0]][dim], iVertices[tri[1]][dim]), iVertices[tri[2]][dim]) + bandDist;
    }

    // Find index ranges of bounding box
    std::array<int, 3> idxMin= {0, 0, 0}, idxMax= {0, 0, 0};
    for (unsigned int dim= 0; dim < 3; dim++) {
      if (nb[dim] > 1) {
        idxMin[dim]= std::max((int)std::round(double(nb[dim] - 1) * (posMin[dim] - start[dim]) * gridBoxSizeInv[dim]), 0);
        idxMax[dim]= std::min((int)std::round(double(nb[dim] - 1) * (posMax[dim] - start[dim]) * gridBoxSizeInv[dim]), nb[dim]-1);
      }
    }

    // Find neighboring voxels that intersect the current triangle
    const Vec::Vec3<double> triFirstVert(iVertices[tri[0]]);
    for (int x= idxMin[0]; x <= idxMax[0]; x++) {
      for (int y= idxMin[1]; y <= idxMax[1]; y++) {
        for (int z= idxMin[2]; z <= idxMax[2]; z++) {
          // Get current point position
          const Vec::Vec3<double> cellPos(startX + double(x) * stepX, startY + double(y) * stepY, startZ + double(z) * stepZ);

          // Quick rejection test based on plane distance
          if (std::abs(faceNormals[k].dot(cellPos-triFirstVert)) > bandDist)
            continue;

          // Get projection on triangle
          Vec::Vec3<double> projPoint;
          int elemType= -1;
          MeshVoxelizer::ComputePointToTriangleProjection(iVertices[tri[0]], iVertices[tri[1]], iVertices[tri[2]], cellPos, projPoint, elemType);
          const double unsignedDist= (cellPos - projPoint).norm();

          // Update distance field if needed
          if (unsignedDist <= bandDist) {
            Vec::Vec3<double> projNormal(0.0, 0.0, 0.0);
            if      (elemType == 0) projNormal= faceNormals[k];
            else if (elemType == 1) projNormal= vertNormals[tri[0]];
            else if (elemType == 2) projNormal= vertNormals[tri[1]];
            else if (elemType == 3) projNormal= vertNormals[tri[2]];
            else if (elemType == 4) projNormal= edgeNormals[{tri[0], tri[1]}];
            else if (elemType == 5) projNormal= edgeNormals[{tri[0], tri[2]}];
            else                    projNormal= edgeNormals[{tri[1], tri[2]}];

            const int xyz= oDist.getFlatIndex(x, y, z);
            if (oDist.at(xyz) > unsignedDist) {
              oDist.at(xyz)= unsignedDist;
              oIsSet.at(xyz)= 1;
              distSign.at(xyz)= ((cellPos - projPoint).dot(projNormal) < 0.0) ? -1 : 1;
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
  // https://www.geometrictools.com/GTE/Mathematics/DistPointTriangle.h
  //         C
  //        /|
  //       / |   
  //   t  /  |
  //     /   |
  //    /____|
  //   A   s  B
  //
  // Proj on face   ABC: oElemType= 0
  // Proj on vertex A:   oElemType= 1
  // Proj on vertex B:   oElemType= 2
  // Proj on vertex C:   oElemType= 3
  // Proj on edge   AB:  oElemType= 4
  // Proj on edge   AC:  oElemType= 5
  // Proj on edge   BC:  oElemType= 6

  // Precomputations
  const Vec::Vec3<double> edge0= iTriangleB - iTriangleA;
  const Vec::Vec3<double> edge1= iTriangleC - iTriangleA;
  const Vec::Vec3<double> diff= iTriangleA - iPoint;
  const double a00= edge0.dot(edge0);
  const double a01= edge0.dot(edge1);
  const double a11= edge1.dot(edge1);
  const double b0= diff.dot(edge0);
  const double b1= diff.dot(edge1);
  const double det= std::max(a00 * a11 - a01 * a01, 0.0);

  double s= a01 * b1 - a11 * b0;
  double t= a01 * b0 - a00 * b1;

  // Barycentric coordinates to find closest point on triangle
  if (s + t <= det) {
    if (s < 0.0) {
      if (t < 0.0) {  // Proj on vertex A:   oElemType= 1
        oElemType= 1;
        if (b0 < 0.0) {
          t= 0.0;
          if (-b0 >= a00) s= 1.0;
          else s= -b0 / a00;
        }
        else {
          s= 0.0;
          if (b1 >= 0.0) t= 0.0;
          else if (-b1 >= a11) t= 1.0;
          else t= -b1 / a11;
        }
      }
      else {  // Proj on edge   AC:  oElemType= 5
        oElemType= 5;
        s= 0.0;
        if (b1 >= 0.0) t= 0.0;
        else if (-b1 >= a11) t= 1.0;
        else t= -b1 / a11;
      }
    }
    else if (t < 0.0) {  // Proj on edge   AB:  oElemType= 4
      oElemType= 4;
      t= 0.0;
      if (b0 >= 0.0) s= 0.0;
      else if (-b0 >= a00) s= 1.0;
      else s= -b0 / a00;
    }
    else {  // Proj on face   ABC: oElemType= 0
      oElemType= 0;
      s/= det;
      t/= det;
    }
  }
  else {
    if (s < 0.0) {  // Proj on vertex C:   oElemType= 3
      oElemType= 3;
      const double tmp0= a01 + b0;
      const double tmp1= a11 + b1;
      if (tmp1 > tmp0) {
        const double numer= tmp1 - tmp0;
        const double denom= a00 - 2.0 * a01 + a11;
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
        else if (b1 >= 0.0) t= 0.0;
        else t= -b1 / a11;
      }
    }
    else if (t < 0.0) {  // Proj on vertex B:   oElemType= 2
      oElemType= 2;
      const double tmp0= a01 + b1;
      const double tmp1= a00 + b0;
      if (tmp1 > tmp0) {
        const double numer= tmp1 - tmp0;
        const double denom= a00 - 2.0 * a01 + a11;
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
        if (tmp1 <= 0.0) s= 1.0;
        else if (b0 >= 0.0) s= 0.0;
        else s= -b0 / a00;
      }
    }
    else {  // Proj on edge   BC:  oElemType= 6
      oElemType= 6;
      const double numer= a11 + b1 - a01 - b0;
      if (numer <= 0.0) {
        s= 0.0;
        t= 1.0;
      }
      else {
        const double denom= a00 - 2.0 * a01 + a11;
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
