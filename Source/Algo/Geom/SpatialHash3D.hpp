#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Vec.hpp"

// Efficient spatial hashing of unbounded space for 3D short-range interactions
// Warning: Hash collisions may occur and lead to double counting. 
// Suggested default settings:
// - The hash table can be the same length as the number of points to add in it
// - The cell size can be the range of the interaction (e.g. 2*radii for collisions)
// 
// Ref document: https://matthias-research.github.io/pages/tenMinutePhysics/11-hashing.pdf
//
// Example usage:
//
// // Allocate and fill the spatial hash with default parameters
// static SpatialHash3D Hash3D;
// const float reach= 2.0 * radius;
// const float reachSqr= reach * reach;
// Hash3D.Fill(Pos, Pos.size(), reach);
// // Prepare the reach offset
// const Vec::Vec3<float> vecOffset(reach);
// // Sweep through the points in parallel
// #pragma omp parallel for
// for (int k0= 0; k0 < N; k0++) {
//   // Sweep through the candidate cells withing reach of the point
//   int minX, minY, minZ, maxX, maxY, maxZ;
//   Hash3D.GetCell(Pos[k0] - vecOffset, minX, minY, minZ);
//   Hash3D.GetCell(Pos[k0] + vecOffset, maxX, maxY, maxZ);
//   // // Mark checked hashes to avoid double checking in case of hash collisions
//   // std::vector<unsigned int> checkedHash;
//   for (int cellX= minX; cellX <= maxX; cellX++) {
//     for (int cellY= minY; cellY <= maxY; cellY++) {
//       for (int cellZ= minZ; cellZ <= maxZ; cellZ++) {
//         // Test all the potential neighbor points in the candidate cells
//         unsigned int hash= Hash3D.GetHash(cellX, cellY, cellZ);
//         // // Skip already checked hashes
//         // if (find(checkedHash.begin(), checkedHash.end(), hash) != checkedHash.end()) continue;
//         // checkedHash.push_back(hash);
//         for (unsigned int idxStreak= Hash3D.beg[hash]; idxStreak <= Hash3D.end[hash]; idxStreak++) {
//           const int k1= Hash3D.idx[idxStreak];
//           if (k0 == k1) continue;
//           // Compute the interaction
//           if ((Pos[k1] - Pos[k0]).normSquared() <= reachSqr) {
//             // Do the thing
//           }
//         }
//       }
//     }
//   }
// }
class SpatialHash3D
{
  private:
  unsigned int nbPoints;
  unsigned int nbHash;
  float cellSize;
  float cellSizeInv;

  public:
  std::vector<unsigned int> beg;
  std::vector<unsigned int> end;
  std::vector<unsigned int> idx;

  SpatialHash3D();
  void Fill(const std::vector<Vec::Vec3<float>>& iArrPos, const unsigned int iNbHash, const float iCellSize);

  inline unsigned int GetFlatIdx(const int iCellX, const int iCellY, const int iCellZ) {
    // Hash function used by https://arxiv.org/pdf/2009.06944
    return ((unsigned int)(iCellX * 73856093) ^
    (unsigned int)(iCellY * 19349663) ^
    (unsigned int)(iCellZ * 83492791)) % nbHash;
  }
  
  inline unsigned int GetFlatIdx(const Vec::Vec3<float>& iPos) {
    int cellX, cellY, cellZ;
    GetCell(iPos, cellX, cellY, cellZ);
    return GetFlatIdx(cellX, cellY, cellZ);
  }
  
  inline void GetCell(const Vec::Vec3<float>& iPos, int& oCellX, int& oCellY, int& oCellZ) {
    oCellX= int(iPos[0] * cellSizeInv);
    oCellY= int(iPos[1] * cellSizeInv);
    oCellZ= int(iPos[2] * cellSizeInv);
  }
};
