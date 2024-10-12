#include "Labelizer.hpp"

// Standard lib
#include <array>
#include <vector>


void Labelizer::LabelConnectedComponents(const int iNbX,
                                         const int iNbY,
                                         const int iNbZ,
                                         const std::vector<int>& iCategoryID,
                                         std::vector<int>& oComponentID) {
  // Check input validity
  if (iNbX <= 0 || iNbY <= 0 || iNbZ <= 0) return;
  if ((int)iCategoryID.size() != iNbX * iNbY * iNbZ) return;

  // Neighbor filter encoding the valid propagation directions
  std::array<std::array<int, 3>, 6> ListNeighborOffsets;
  ListNeighborOffsets[0]= {-1, +0, +0};
  ListNeighborOffsets[1]= {+1, +0, +0};
  ListNeighborOffsets[2]= {+0, -1, +0};
  ListNeighborOffsets[3]= {+0, +1, +0};
  ListNeighborOffsets[4]= {+0, +0, -1};
  ListNeighborOffsets[5]= {+0, +0, +1};

  // Initialize label field, queue and current label
  oComponentID.resize(iNbX * iNbY * iNbZ);
  std::fill(oComponentID.begin(), oComponentID.end(), -1);
  std::vector<std::array<int, 3>> Q;
  int currentLabel= 0;

  // Sweep through the field
  for (int xCheck= 0; xCheck < iNbX; xCheck++) {
    for (int yCheck= 0; yCheck < iNbY; yCheck++) {
      for (int zCheck= 0; zCheck < iNbZ; zCheck++) {
        // Find unlabeled voxel
        const int xyzCheck= xCheck * iNbY * iNbZ + yCheck * iNbZ + zCheck;
        if (oComponentID[xyzCheck] < 0) {
          Q.push_back({xCheck, yCheck, zCheck});

          // Propagate the label through the entire connected component
          while (!Q.empty()) {
            const int xCurr= Q[Q.size() - 1][0];
            const int yCurr= Q[Q.size() - 1][1];
            const int zCurr= Q[Q.size() - 1][2];
            Q.pop_back();
            const int xyzCurr= xCurr * iNbY * iNbZ + yCurr * iNbZ + zCurr;
            oComponentID[xyzCurr]= currentLabel;

            // Look for unlabeled voxels in the immediate neighbors
            for (std::array<int, 3> neighborOffset : ListNeighborOffsets) {
              const int xProp= xCurr + neighborOffset[0];
              const int yProp= yCurr + neighborOffset[1];
              const int zProp= zCurr + neighborOffset[2];
              if (xProp >= 0 && yProp >= 0 && zProp >= 0 && xProp < iNbX && yProp < iNbY && zProp < iNbZ) {
                const int xyzProp= xProp * iNbY * iNbZ + yProp * iNbZ + zProp;
                if (oComponentID[xyzProp] < 0 && iCategoryID[xyzCurr] == iCategoryID[xyzProp]) {
                  Q.push_back({xProp, yProp, zProp});
                }
              }
            }
          }

          // Increase the current label before starting to label the next component
          currentLabel++;
        }
      }
    }
  }
}
