#pragma once

// Standard lib
#include <vector>


class Labelizer
{
  public:
  // Compute a unique label to each connected component in 2D sharing a common category ID
  // Uses 4-neighborhood connectivity, generated labels are positive and start at 0
  // Input : 2D field with different integer value for each category of voxel
  // Output: unique integer for each connected component
  // Returns the number of connected components
  static int LabelConnectedComponents(const int iNbX,
                                      const int iNbY,
                                      const std::vector<int>& iCategoryID,
                                      std::vector<int>& oComponentID);

  // Compute a unique label to each connected component in 3D sharing a common category ID
  // Uses 6-neighborhood connectivity, generated labels are positive and start at 0
  // Input : 3D field with different integer value for each category of voxel
  // Output: unique integer for each connected component
  // Returns the number of connected components
  static int LabelConnectedComponents(const int iNbX,
                                      const int iNbY,
                                      const int iNbZ,
                                      const std::vector<int>& iCategoryID,
                                      std::vector<int>& oComponentID);
};
