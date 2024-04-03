#include "Sketch.hpp"

// Standard lib
#include <cmath>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


void Sketch::PolylineSubdivideAndSmooth(
    bool const iIsOpenPolyline,
    int const iNbNewNodesPerSegment,
    int const iNbSmoothingSteps,
    std::vector<std::array<double, 3>>& ioPolyline) {
  // Create a new list of vertices
  std::vector<Vec::Vec3<double>> vertices;

  // Subdivide each segment with linear interpolation
  if (iIsOpenPolyline) {
    for (int k0= 0; k0 < int(ioPolyline.size()) - 1; k0++) {
      vertices.push_back(ioPolyline[k0]);
      for (int k1= 0; k1 < iNbNewNodesPerSegment; k1++) {
        double ratio= double(k1 + 1) / double(iNbNewNodesPerSegment + 1);
        vertices.push_back((1.0 - ratio) * Vec::Vec3<double>(ioPolyline[k0]) + ratio * Vec::Vec3<double>(ioPolyline[k0 + 1]));
      }
    }
    vertices.push_back(ioPolyline[int(ioPolyline.size()) - 1]);
  }
  else {
    for (int k0= 0; k0 < int(ioPolyline.size()); k0++) {
      vertices.push_back(ioPolyline[k0]);
      for (int k1= 0; k1 < iNbNewNodesPerSegment; k1++) {
        double ratio= double(k1 + 1) / double(iNbNewNodesPerSegment + 1);
        vertices.push_back((1.0 - ratio) * Vec::Vec3<double>(ioPolyline[k0]) + ratio * Vec::Vec3<double>(ioPolyline[(k0 + 1) % int(ioPolyline.size())]));
      }
    }
  }

  // Apply Laplacian smoothing
  for (int idxSmoothStep= 0; idxSmoothStep < iNbSmoothingSteps; idxSmoothStep++) {
    std::vector<Vec::Vec3<double>> verticesOld= vertices;
    int const N= int(vertices.size());
    if (iIsOpenPolyline)
      for (int k= 1; k < int(vertices.size()) - 1; k++)
        vertices[k]= (verticesOld[k - 1] + verticesOld[k] + verticesOld[k + 1]) / 3.0;
    else
      for (int k= 0; k < N; k++)
        vertices[k]= (verticesOld[(k - 1 + N) % N] + verticesOld[k] + verticesOld[(k + 1) % N]) / 3.0;
  }

  // Output the new polyline
  ioPolyline.clear();
  for (int k0= 0; k0 < int(vertices.size()); k0++)
    ioPolyline.push_back({vertices[k0][0], vertices[k0][1], vertices[k0][2]});
}
