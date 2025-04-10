#include "SpatialHash3D.hpp"

// #define SPATIALHASH3D_DEBUG_PRINT_ARRAY

// Standard lib
#include <cassert>
#include <vector>
#ifdef SPATIALHASH3D_DEBUG_PRINT_ARRAY
#include <cstdio>
#endif

// Algo headers
#include "Type/Vec.hpp"


SpatialHash3D::SpatialHash3D() {
  cellSize= cellSizeInv= 1.0f;
  nbHash= 1;
  nbPoints= 1;
  beg.resize(1);
  end.resize(1);
  idx.resize(1);
}


void SpatialHash3D::Fill(const std::vector<Vec::Vec3<float>>& iArrPos, const unsigned int iNbHash, const float iCellSize) {
  // Check the input parameters
  assert(!iArrPos.empty() && iNbHash > 0  && iCellSize > 0.0f);

  cellSize= iCellSize;
  cellSizeInv= 1.0f / iCellSize;
  nbHash= iNbHash;
  nbPoints= iArrPos.size();
  beg.resize(nbHash);
  end.resize(nbHash);
  idx.resize(nbPoints);

  // Compute the counts of points for each hash
  std::fill(end.begin(), end.end(), 0);        // Initialize the counts to 0
  for (unsigned int k= 0; k < nbPoints; k++) { // Sweep through the points
    end[GetFlatIdx(iArrPos[k])]++;                // Increase the count for the hash of the current point
  }
  #ifdef SPATIALHASH3D_DEBUG_PRINT_ARRAY
  printf("siz: "); for (unsigned int k : end) printf("%d ", k); printf("\n");
  #endif

  // Compute the cumulative sum of the counts
  // Compute the index of the last point for each hash streak
  unsigned int cumSum= 0;                              // Initialize the cumulative sum value to 0
  for (unsigned int hash= 0; hash < nbHash; hash++) {  // Sweep through the potential hash values
    cumSum+= end[hash];                                // Increase the cumulative sum value
    beg[hash]= cumSum;                                 // Set beg as a cumulative sum of the counts
    end[hash]= (beg[hash] > 0) * (beg[hash] - 1);      // Set end as the index of the last point 
  }

  // Compute the index of the first point for each hash streak
  // Insert each point index in the correct position
  for (unsigned int k= 0; k < nbPoints; k++) {    // Sweep through the points
    const unsigned int hash= GetFlatIdx(iArrPos[k]); // Get the hash
    beg[hash]--;                                  // Decrease the cumulative count down to the index of the beginning
    idx[beg[hash]]= k;                            // Insert the point index in its position
  }
  #ifdef SPATIALHASH3D_DEBUG_PRINT_ARRAY
  printf("beg: "); for (unsigned int k : beg) printf("%d ", k); printf("\n");
  printf("end: "); for (unsigned int k : end) printf("%d ", k); printf("\n");
  printf("idx: "); for (unsigned int k : idx) printf("%d ", k); printf("\n");
  #endif
}
