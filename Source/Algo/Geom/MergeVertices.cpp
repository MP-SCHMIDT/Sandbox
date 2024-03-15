#include "MergeVertices.hpp"

// Standard lib
#include <array>
#include <vector>


void MergeVertices::QuadraticMerge(
    double const iTolerance,
    std::vector<std::array<double, 3>>& ioVertices,
    std::vector<std::array<int, 3>>& ioTriangles) {
  // Initialize the label vector with original indices
  std::vector<int> label(ioVertices.size());
  for (int vertIdx= 0; vertIdx < (int)ioVertices.size(); vertIdx++)
    label[vertIdx]= vertIdx;

  // Spread the smallest label through neighbors when distance is smaller than iTolerance
  bool hasSpread= true;
  while (hasSpread) {
    hasSpread= false;
    for (int k0= 0; k0 < (int)ioVertices.size(); k0++) {
      for (int k1= 0; k1 < (int)ioVertices.size(); k1++) {
        const double distSqr= (ioVertices[k0][0] - ioVertices[k1][0]) * (ioVertices[k0][0] - ioVertices[k1][0]) +
                              (ioVertices[k0][1] - ioVertices[k1][1]) * (ioVertices[k0][1] - ioVertices[k1][1]) +
                              (ioVertices[k0][2] - ioVertices[k1][2]) * (ioVertices[k0][2] - ioVertices[k1][2]);
        if (k0 != k1 && distSqr <= iTolerance * iTolerance) {
          if (label[k0] > label[k1]) {
            label[k0]= label[k1];
            hasSpread= true;
          }
        }
      }
    }
  }

  // Create the compressed list of vertices and the associated index mappings
  std::vector<std::array<double, 3>> oldVertices= ioVertices;
  std::vector<int> idxMap(ioVertices.size());
  ioVertices.clear();
  for (int idxVert= 0; idxVert < (int)oldVertices.size(); idxVert++) {
    if (idxVert == label[idxVert]) {
      label[idxVert]= (int)ioVertices.size();
      ioVertices.push_back(oldVertices[idxVert]);
      idxMap[idxVert]= label[idxVert];
    }
    else {
      idxMap[idxVert]= label[label[idxVert]];
    }
  }

  // Create the new triangles using the merged index mappings
  std::vector<std::array<int, 3>> oldTriangles= ioTriangles;
  ioTriangles.clear();
  for (int idxTri= 0; idxTri < (int)oldTriangles.size(); idxTri++) {
    if (idxMap[oldTriangles[idxTri][0]] == idxMap[oldTriangles[idxTri][1]]) continue;
    if (idxMap[oldTriangles[idxTri][0]] == idxMap[oldTriangles[idxTri][2]]) continue;
    if (idxMap[oldTriangles[idxTri][1]] == idxMap[oldTriangles[idxTri][2]]) continue;
    ioTriangles.push_back({idxMap[oldTriangles[idxTri][0]], idxMap[oldTriangles[idxTri][1]], idxMap[oldTriangles[idxTri][2]]});
  }
}
