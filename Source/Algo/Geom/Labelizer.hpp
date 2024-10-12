#pragma once

// Standard lib
#include <vector>


class Labelizer
{
  public:
  // Compute a unique label to each connected component sharing a common category ID
  // Uses 6-neighborhood connectivity, generated labels are positive and start at 0
  // Input : field with different integer value for each category of voxel
  // Output: unique integer for each connected component
  static void LabelConnectedComponents(const int iNbX,
                                       const int iNbY,
                                       const int iNbZ,
                                       const std::vector<int>& iCategoryID,
                                       std::vector<int>& oComponentID);
};
