#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Type/Field.hpp"
#include "Type/Vec.hpp"


class MeshVoxelizer
{
  public:
  // Computes a narrow band signed distance fields to a triangular meshes
  // iVertices    : Mesh vertex positions
  // iTriangles   : Mesh triangles index triplets
  // iBoxMin      : Bounding box of the distance field
  // iBoxMax      : Bounding box of the distance field
  // iNbX         : Resolution of the distance field
  // iNbY         : Resolution of the distance field
  // iNbZ         : Resolution of the distance field
  // iIsCentered  : Flag for center-valued or corner-valued voxels
  // oIsSet       : Flag field defining if the voxel is in the narrow band
  // oDistSigned  : Distance value
  static void ComputeNarrowBandSignedDistanceField(const std::vector<std::array<double, 3>>& iVertices,
                                                   const std::vector<std::array<int, 3>>& iTriangles,
                                                   const std::array<double, 3>& iBoxMin,
                                                   const std::array<double, 3>& iBoxMax,
                                                   const int iNbX,
                                                   const int iNbY,
                                                   const int iNbZ,
                                                   const bool iIsCentered,
                                                   Field::Field3<char>& oIsSet,
                                                   Field::Field3<double>& oDist);


  // Utility function used in the computation of narrow band signed distance fields to triangular meshes
  // Projects a given 3D point to the closest point of a 3D triangle
  // The integer type informs if the projected point landed on the face (0), on an edge (1, 2, 3) or on a corner (4, 5, 6)
  static void ComputePointToTriangleProjection(const Vec::Vec3<double>& iTriangleA,
                                               const Vec::Vec3<double>& iTriangleB,
                                               const Vec::Vec3<double>& iTriangleC,
                                               const Vec::Vec3<double>& iPoint,
                                               Vec::Vec3<double>& oProjPoint,
                                               int& oElemType);
};
