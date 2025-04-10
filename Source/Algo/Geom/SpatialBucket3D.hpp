#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Type/Vec.hpp"

// Efficient spatial partition in fixed size buckets for 3D short-range interactions (e.g. collisions)
// Suggested default settings:
// - The bucket size can be the range of the interaction (e.g. 2*radii for collisions)
// - The bucket capacity can be set to e.g. 10 and increased in case of overflow
// - If the max buckets count is specified, the bucket size is automatically ajusted keep the number of buckets in range
//
// Example usage:
// Bucket3D.SetDim(boxMin, boxMax, reach, bucketCapacity, maxBucketCount);
// BucketOverflown= Bucket3D.Fill(Pos);
// // Sweep through the points
// const Vec::Vec3<float> vecOffset(maxForceLawRange);
// #pragma omp parallel for
// for (int k0= 0; k0 < (int)Pos.size(); k0++) {
//   // Sweep through the candidate cells withing reach of the point
//   int minX, minY, minZ, maxX, maxY, maxZ;
//   Bucket3D.GetCell(Pos[k0] - vecOffset, minX, minY, minZ);
//   Bucket3D.GetCell(Pos[k0] + vecOffset, maxX, maxY, maxZ);
//   for (int cellX= std::max(0, minX); cellX <= std::min(maxX, Bucket3D.nX - 1); cellX++) {
//     for (int cellY= std::max(0, minY); cellY <= std::min(maxY, Bucket3D.nY - 1); cellY++) {
//       for (int cellZ= std::max(0, minZ); cellZ <= std::min(maxZ, Bucket3D.nZ - 1); cellZ++) {
//         // Test all the potential neighbor points in the candidate cells
//         for (int k1 : Bucket3D.buckets[Bucket3D.GetFlatIdx(cellX, cellY, cellZ)]) {
//           if (k0 == k1) continue;
class SpatialBucket3D
{
  private:
  float stepXinv;
  float stepYinv;
  float stepZinv;
  int nYZ;
  
  public:
  int nX;
  int nY;
  int nZ;
  int nXYZ;
  int bucketCapacity;
  Vec::Vec3<float> boxMin;
  Vec::Vec3<float> boxMax;
  std::vector<std::vector<unsigned int>> buckets;

  SpatialBucket3D();
  void SetDim(const Vec::Vec3<float>& iBoxMin, const Vec::Vec3<float>& iBoxMax, const float iMinBucketSize, const int iBucketCapacity, const int iMaxNbBuckets= 0);
  bool Fill(const std::vector<Vec::Vec3<float>>& iArrPos);
  
  inline void GetCell(const Vec::Vec3<float>& iPos, int& oIdxX, int& oIdxY, int& oIdxZ) {
    oIdxX= int(std::floor((iPos[0] - boxMin[0]) * stepXinv));
    oIdxY= int(std::floor((iPos[1] - boxMin[1]) * stepYinv));
    oIdxZ= int(std::floor((iPos[2] - boxMin[2]) * stepZinv));
  }

  inline int GetFlatIdx(const int iIdxX, const int iIdxY, const int iIdxZ) {
    return iIdxX * nYZ + iIdxY*nZ + iIdxZ;
  }

};
