#pragma once

// Standard lib
#include <array>
#include <cmath>


namespace BoxGrid {
  inline void GetVoxelSizes(
      int const iNbX,
      int const iNbY,
      int const iNbZ,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      bool const iCentered,
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
      int const iNbX,
      int const iNbY,
      int const iNbZ,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      bool const iCentered,
      double& oVoxSizeX,
      double& oVoxSizeY,
      double& oVoxSizeZ,
      double& oVoxSizeDiag) {
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
    oVoxSizeDiag= std::sqrt(oVoxSizeX * oVoxSizeX + oVoxSizeY * oVoxSizeY + oVoxSizeZ * oVoxSizeZ);
  }


  inline void GetVoxelStart(
      std::array<double, 3> const& iBBoxMin,
      double const iVoxSizeX,
      double const iVoxSizeY,
      double const iVoxSizeZ,
      bool const iCentered,
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
