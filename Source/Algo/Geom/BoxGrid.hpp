#pragma once

// Standard lib
#include <array>
#include <cmath>


namespace BoxGrid {
  inline void GetVoxelSizes(
      const int iNbX,
      const int iNbY,
      const int iNbZ,
      const std::array<double, 3>& iBBoxMin,
      const std::array<double, 3>& iBBoxMax,
      const bool iCentered,
      double& oVoxSizeX,
      double& oVoxSizeY,
      double& oVoxSizeZ) {
    if (iCentered) {
      oVoxSizeX= (iBBoxMax[0] - iBBoxMin[0]) / double(iNbX);
      oVoxSizeY= (iBBoxMax[1] - iBBoxMin[1]) / double(iNbY);
      oVoxSizeZ= (iBBoxMax[2] - iBBoxMin[2]) / double(iNbZ);
    }
    else {
      oVoxSizeX= (iBBoxMax[0] - iBBoxMin[0]) / double(iNbX - 1);
      oVoxSizeY= (iBBoxMax[1] - iBBoxMin[1]) / double(iNbY - 1);
      oVoxSizeZ= (iBBoxMax[2] - iBBoxMin[2]) / double(iNbZ - 1);
    }
  }

  
  inline void GetVoxelSizes(
      const int iNbX,
      const int iNbY,
      const int iNbZ,
      const std::array<double, 3>& iBBoxMin,
      const std::array<double, 3>& iBBoxMax,
      const bool iCentered,
      double& oVoxSizeX,
      double& oVoxSizeY,
      double& oVoxSizeZ,
      double& oVoxSizeDiag) {
    GetVoxelSizes(iNbX, iNbY, iNbZ, iBBoxMin, iBBoxMax, iCentered, oVoxSizeX, oVoxSizeY, oVoxSizeZ);
    oVoxSizeDiag= std::sqrt(oVoxSizeX * oVoxSizeX + oVoxSizeY * oVoxSizeY + oVoxSizeZ * oVoxSizeZ);
  }


  inline void GetVoxelStart(
      const std::array<double, 3>& iBBoxMin,
      const double iVoxSizeX,
      const double iVoxSizeY,
      const double iVoxSizeZ,
      const bool iCentered,
      double& oStartX,
      double& oStartY,
      double& oStartZ) {
    if (iCentered) {
      oStartX= 0.5 * iVoxSizeX + iBBoxMin[0];
      oStartY= 0.5 * iVoxSizeY + iBBoxMin[1];
      oStartZ= 0.5 * iVoxSizeZ + iBBoxMin[2];
    }
    else {
      oStartX= iBBoxMin[0];
      oStartY= iBBoxMin[1];
      oStartZ= iBBoxMin[2];
    }
  }

}  // namespace BoxGrid
