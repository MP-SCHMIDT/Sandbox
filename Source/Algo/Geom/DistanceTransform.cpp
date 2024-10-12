
#include "DistanceTransform.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

// Algo headers
#include "Type/Field.hpp"


void DistanceTransform::ApplyDistanceTransform(
    const int iMaxIterConvergence,
    const float iSizeX,
    const float iSizeY,
    const float iSizeZ,
    const Field::Field3<char> &iVoxType,
    Field::Field3<float> &ioDist,
    const int iVerbose) {
  // Check inputs
  if (iMaxIterConvergence <= 0 ||
      iVoxType.nXYZ != ioDist.nXYZ ||
      iVoxType.nXYZ <= 0) {
    printf("[ERROR]: Invalid inputs for DistanceTransform::ApplyDistanceTransform\n");
    return;
  }

  // Initialize variables
  const int nX= iVoxType.nX;
  const int nY= iVoxType.nY;
  const int nZ= iVoxType.nZ;
  Field::Field3<char> isSet(nX, nY, nZ, 0);
  for (int xyz= 0; xyz < nX * nY * nZ; xyz++)
    if (iVoxType.at(xyz) > 0)
      isSet.at(xyz)= 1;

  // Build 3D neighbor filter encoding the distance transform for the forward and backward passes
  std::vector<std::array<int, 3>> Neighbor;
  std::vector<float> NeighborDist;
  for (int dx= -1; dx <= 1; dx++) {
    for (int dy= -1; dy <= 1; dy++) {
      for (int dz= -1; dz <= 1; dz++) {
        if (nX == 1 && dx != 0) continue;
        if (nY == 1 && dy != 0) continue;
        if (nZ == 1 && dz != 0) continue;
        if ((dx < 0) || (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz < 0)) {
          Neighbor.push_back({dx, dy, dz});
          NeighborDist.push_back(std::sqrt(float(dx * dx * iSizeX * iSizeX + dy * dy * iSizeY * iSizeY + dz * dz * iSizeZ * iSizeZ)));
        }
      }
    }
  }

  // Iterate the distance transform util convergence or max iterations reached
  if (iVerbose > 0) printf("Iter: ");
  bool isConverged= false;
  for (int iter= 0; iter < iMaxIterConvergence && !isConverged; iter++) {
    if (iVerbose > 0) printf("%d ", iter);
    isConverged= true;
    // Execute the forward then backward pass of the distance transform
    for (int idxPass= 0; idxPass < 2; idxPass++) {
      // Set the loop settings for the current pass
      const int xBeg= (idxPass == 0) ? 0 : nX - 1;
      const int yBeg= (idxPass == 0) ? 0 : nY - 1;
      const int zBeg= (idxPass == 0) ? 0 : nZ - 1;
      const int xEnd= (idxPass == 0) ? nX : -1;
      const int yEnd= (idxPass == 0) ? nY : -1;
      const int zEnd= (idxPass == 0) ? nZ : -1;
      const int xInc= (idxPass == 0) ? 1 : -1;
      const int yInc= (idxPass == 0) ? 1 : -1;
      const int zInc= (idxPass == 0) ? 1 : -1;
      // Sweep through the field
      for (int x= xBeg; x != xEnd; x+= xInc) {
        for (int y= yBeg; y != yEnd; y+= yInc) {
          for (int z= zBeg; z != zEnd; z+= zInc) {
            if (iVoxType.at(x, y, z) != 0) continue;
            for (int k= 0; k < int(Neighbor.size()); k++) {
              const int xShift= (idxPass == 0) ? (x + Neighbor[k][0]) : (x - Neighbor[k][0]);
              const int yShift= (idxPass == 0) ? (y + Neighbor[k][1]) : (y - Neighbor[k][1]);
              const int zShift= (idxPass == 0) ? (z + Neighbor[k][2]) : (z - Neighbor[k][2]);
              if (xShift < 0 || xShift > nX - 1) continue;
              if (yShift < 0 || yShift > nY - 1) continue;
              if (zShift < 0 || zShift > nZ - 1) continue;
              if (iVoxType.at(xShift, yShift, zShift) < 0) continue;
              if (!isSet.at(xShift, yShift, zShift)) continue;
              if (!isSet.at(x, y, z) || std::abs(ioDist.at(xShift, yShift, zShift)) + NeighborDist[k] < std::abs(ioDist.at(x, y, z))) {
                isSet.at(x, y, z)= 1;
                isConverged= false;
                if (ioDist.at(xShift, yShift, zShift) < 0.0f)
                  ioDist.at(x, y, z)= ioDist.at(xShift, yShift, zShift) - NeighborDist[k];
                else
                  ioDist.at(x, y, z)= ioDist.at(xShift, yShift, zShift) + NeighborDist[k];
              }
            }
          }
        }
      }
    }
  }
  if (iVerbose > 0) printf("Converged=%c\n", isConverged ? 'Y' : 'N');
}
