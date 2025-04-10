#include "NBodyGravDynam.hpp"


// Standard lib
#include <cmath>
#include <numeric>
#include <vector>

// Algo headers
#include "Type/Vec.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// TODO try parallelizing tree construction https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
// TODO try compacting the tree by removing empty cells ? Tricky to properly re-wire child/next indices. Might only be possible if tree is in perfect depth-first order before sorting ?
void NBodyGravDynam::BuildTree(const std::vector<Vec::Vec3<float>>& iPos) {
  // Get UI parameters
  const bool infiniteBox= D.UI[TreeInfiniteBox_].I() == 1;
  const bool torusPos= D.UI[DomainTorusPos__].I() == 1;
  const unsigned int treeMaxDepth= std::max(D.UI[TreeMaxDepth____].I(), 0);

  // Define the root cell center and size
  Vec::Vec3<float> rootCenter(0.5f, 0.5f, 0.5f);
  float rootSize= 1.0f;
  if (infiniteBox) {
    Vec::Vec3<float> PosBoxMin= iPos[0], PosBoxMax= iPos[0];
    for (unsigned int k0= 1; k0 < N; k0++) {
      PosBoxMin= PosBoxMin.cwiseMin(iPos[k0]);
      PosBoxMax= PosBoxMax.cwiseMax(iPos[k0]);
    }
    rootCenter= 0.5f * (PosBoxMin + PosBoxMax);
    rootSize= (PosBoxMax - PosBoxMin).maxCoeff();
  }

  // Initialize the tree with the root cell
  Tree.clear();
  Tree.push_back({rootSize, rootCenter, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, 0});

  // Sequentially add the particles in the tree, adding cells as required
  for (unsigned int k0= 0; k0 < N; k0++) {
    if (!infiniteBox && !torusPos && (iPos[k0][0] < 0.0f || iPos[k0][0] > 1.0f ||
                                      iPos[k0][1] < 0.0f || iPos[k0][1] > 1.0f ||
                                      iPos[k0][2] < 0.0f || iPos[k0][2] > 1.0f)) continue;
    // Start the tree search at the root cell
    unsigned int idxCell= 0, curDepth= 0;
    while (true) {
      // Case for empty leaf cell or a non-empty max-depth leaf cell
      if (Tree[idxCell].Count == 0 || curDepth == treeMaxDepth) {
        Tree[idxCell].AvgPos+= iPos[k0];
        Tree[idxCell].AvgVel+= Vel[k0];
        Tree[idxCell].Count++;
        break;
      }
      // Case for leaf cell with no children (by induction, not at max depth and must contain exactly one body)
      if (Tree[idxCell].Child == 0) {
        // Add the 8 children nodes at the end of the array
        const float childSize= 0.5f * Tree[idxCell].Size;
        const float childShift= 0.5f * childSize;
        const unsigned int idxFirstChild= Tree.size();
        const Vec::Vec3<float> iCenter= Tree[idxCell].Center;
        const float c0Nega= iCenter[0] - childShift, c1Nega= iCenter[1] - childShift, c2Nega= iCenter[2] - childShift;
        const float c0Posi= iCenter[0] + childShift, c1Posi= iCenter[1] + childShift, c2Posi= iCenter[2] + childShift;
        Tree.insert(Tree.end(), {
          {childSize, {c0Nega, c1Nega, c2Nega}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+1},
          {childSize, {c0Nega, c1Nega, c2Posi}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+2},
          {childSize, {c0Nega, c1Posi, c2Nega}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+3},
          {childSize, {c0Nega, c1Posi, c2Posi}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+4},
          {childSize, {c0Posi, c1Nega, c2Nega}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+5},
          {childSize, {c0Posi, c1Nega, c2Posi}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+6},
          {childSize, {c0Posi, c1Posi, c2Nega}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, idxFirstChild+7},
          {childSize, {c0Posi, c1Posi, c2Posi}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, 0, 0, Tree[idxCell].Next}
        });
        Tree[idxCell].Child= idxFirstChild;

        // Reinsert the existing body from the current cell into the correct sub cell
        const unsigned int idxChild= PickSubCell(Tree[idxCell].Center, Tree[idxCell].AvgPos, Tree[idxCell].Child);
        Tree[idxChild].AvgPos+= Tree[idxCell].AvgPos;
        Tree[idxChild].AvgVel+= Tree[idxCell].AvgVel;
        Tree[idxChild].Count++;
      }

      // Add the body to the cell and keep going down the tree
      Tree[idxCell].AvgPos+= iPos[k0];
      Tree[idxCell].AvgVel+= Vel[k0];
      Tree[idxCell].Count++;
  
      // Get the child cell index for the new body and set as reference for the next pass in the loop
      idxCell= PickSubCell(Tree[idxCell].Center, iPos[k0], Tree[idxCell].Child);
      curDepth++;
    }
  }

  // Compute the actual center of mass for each tree cell
  for (unsigned int idxCell= 0; idxCell < Tree.size(); idxCell++) {
    if (Tree[idxCell].Count > 0) {
      Tree[idxCell].AvgPos= Tree[idxCell].AvgPos/float(Tree[idxCell].Count);
      Tree[idxCell].AvgVel= Tree[idxCell].AvgVel/float(Tree[idxCell].Count);
    }
  }
}
