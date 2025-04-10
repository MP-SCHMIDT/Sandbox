#include "NBodyGravDynam.hpp"


// Standard lib
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// Algo headers
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Sort the bodies in Mordon order (a.k.a. Z-order curve)
void NBodyGravDynam::MortonSortBodies() {
  // Compute the bounding box
  Vec::Vec3<float> PosBoxMin= Pos[0], PosBoxMax= Pos[0];
  for (unsigned int k0= 1; k0 < N; k0++) {
    PosBoxMin= PosBoxMin.cwiseMin(Pos[k0]);
    PosBoxMax= PosBoxMax.cwiseMax(Pos[k0]);
  }

  // Compute the Morton index for each body
  std::vector<unsigned int> mortonIndex(Pos.size());
  for (unsigned int k0= 0; k0 < N; k0++) {
    const Vec::Vec3<float> unitPos= (Pos[k0]-PosBoxMin).cwiseDiv(PosBoxMax-PosBoxMin);
    mortonIndex[k0]= MortonGetIndex(unitPos[0], unitPos[1], unitPos[2]);
  }

  // Sort the bodies positions and velocities
  std::vector<unsigned int> idx(Pos.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](const unsigned int A, const unsigned int B) -> bool { return mortonIndex[A] < mortonIndex[B]; });
  const std::vector<Vec::Vec3<float>> PosOld= Pos;
  const std::vector<Vec::Vec3<float>> VelOld= Vel;
  for (unsigned int k0= 0; k0 < N; k0++) {
    Pos[k0]= PosOld[idx[k0]];
    Vel[k0]= VelOld[idx[k0]];
  }
}


// Expands a 10-bit integer into 30 bits by inserting 2 zeros after each bit.
// https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
unsigned int NBodyGravDynam::MortonExpandBits(unsigned int v) {
    v= (v * 0x00010001u) & 0xFF0000FFu;
    v= (v * 0x00000101u) & 0x0F00F00Fu;
    v= (v * 0x00000011u) & 0xC30C30C3u;
    v= (v * 0x00000005u) & 0x49249249u;
    return v;
}


// Calculates a 30-bit Morton code for the given 3D point located within the unit cube [0,1].
// https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
unsigned int NBodyGravDynam::MortonGetIndex(float x, float y, float z) {
    x= std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
    y= std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
    z= std::min(std::max(z * 1024.0f, 0.0f), 1023.0f);
    unsigned int xx= MortonExpandBits((unsigned int)x);
    unsigned int yy= MortonExpandBits((unsigned int)y);
    unsigned int zz= MortonExpandBits((unsigned int)z);
    return xx * 4 + yy * 2 + zz;
}
