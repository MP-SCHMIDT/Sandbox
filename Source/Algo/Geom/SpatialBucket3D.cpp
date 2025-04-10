#include "SpatialBucket3D.hpp"

// Standard lib
#include <cassert>
#include <vector>
#include <cstdio>

// Algo headers
#include "Type/Vec.hpp"


SpatialBucket3D::SpatialBucket3D() {
  buckets= std::vector<std::vector<unsigned int>>(1, std::vector<unsigned int>(1, 0));
}


void SpatialBucket3D::SetDim(const Vec::Vec3<float>& iBoxMin,
                             const Vec::Vec3<float>& iBoxMax,
                             const float iMinBucketSize,
                             const int iBucketCapacity,
                             const int iMaxNbBuckets) {
  // Check the input parameters
  assert(iMinBucketSize > 0.0f && iBucketCapacity > 0 && (iBoxMax - iBoxMin).minCoeff() >= 0.0f);

  // Compute the grid dimensions
  const float boxDX= iBoxMax[0] - iBoxMin[0];
  const float boxDY= iBoxMax[1] - iBoxMin[1];
  const float boxDZ= iBoxMax[2] - iBoxMin[2];
  float bucketSize= iMinBucketSize;
  nX= std::max(1, int(std::ceil(boxDX / bucketSize)));
  nY= std::max(1, int(std::ceil(boxDY / bucketSize)));
  nZ= std::max(1, int(std::ceil(boxDZ / bucketSize)));

  // Optimize the bucket size to get the smallest size above the threshold and yielding less buckets than the upper limit
  if (iMaxNbBuckets > 0) {
    float bucketSizeTest= bucketSize;
    bucketSize= std::max(std::max(boxDX, boxDY), boxDZ);
    for (int iter= 0; iter < 20; iter++) {
      const int nXtest= std::max(1, int(std::ceil(boxDX / bucketSizeTest)));
      const int nYtest= std::max(1, int(std::ceil(boxDY / bucketSizeTest)));
      const int nZtest= std::max(1, int(std::ceil(boxDZ / bucketSizeTest)));
      if (nXtest * nYtest * nZtest > iMaxNbBuckets || bucketSizeTest < iMinBucketSize) {
        // Increase the bucket size to have fewer of them or get above the provided threshold
        bucketSizeTest*= 1.5f;
      }
      else {
        // Save the current 
        if (bucketSizeTest < bucketSize) {
          nX= nXtest;
          nY= nYtest;
          nZ= nZtest;
          bucketSize= bucketSizeTest;
        }
        bucketSizeTest*= 0.5f;
      }
    }
  }

  // Set the final dimensions
  nYZ= nY * nZ;
  nXYZ= nX * nYZ;
  bucketCapacity= iBucketCapacity;
  const Vec::Vec3<float> boxSize(float(nX) * bucketSize, float(nY) * bucketSize, float(nZ) * bucketSize);
  const Vec::Vec3<float> boxCenter(0.5f * (iBoxMin + iBoxMax));
  boxMin= boxCenter - 0.5 * boxSize;
  boxMax= boxCenter + 0.5 * boxSize;
  stepXinv= 1.0f / (bucketSize * 1.001f); // Tiny scaling to transparenlty handle points exactly on the boxMax face of the domain
  stepYinv= 1.0f / (bucketSize * 1.001f); // Tiny scaling to transparenlty handle points exactly on the boxMax face of the domain
  stepZinv= 1.0f / (bucketSize * 1.001f); // Tiny scaling to transparenlty handle points exactly on the boxMax face of the domain
}


bool SpatialBucket3D::Fill(const std::vector<Vec::Vec3<float>>& iArrPos) {
  // Check the input parameters
  assert(!iArrPos.empty());

  // Allocate the buckets if necessary
  if ((int)buckets.size() != nXYZ || (int)buckets[0].capacity() != bucketCapacity)
    buckets= std::vector<std::vector<unsigned int>>(nXYZ, std::vector<unsigned int>(bucketCapacity));
  for (int idxBucket= 0; idxBucket < nXYZ; idxBucket++)
    buckets.at(idxBucket).clear();
  
  // Fill the spatial partition
  bool hasOverflown= false;
  for (unsigned int k= 0; k < iArrPos.size(); k++) {
    int idxX, idxY, idxZ;
    GetCell(iArrPos[k], idxX, idxY, idxZ);
    if (idxX >= 0 && idxX < nX && idxY >= 0 && idxY < nY && idxZ >= 0 && idxZ < nZ) {
      const int idxXYZ= GetFlatIdx(idxX, idxY, idxZ);
      if ((int)buckets[idxXYZ].size() < bucketCapacity)
        buckets[idxXYZ].push_back(k);
      else
        hasOverflown= true;
    }
  }

  return hasOverflown;
}
