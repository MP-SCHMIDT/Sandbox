#include "Bresenham.hpp"

// Standard lib
#include <array>
#include <cmath>
#include <vector>


std::vector<std::array<int, 2>> Bresenham::Line2D(int x0, int y0, int x1, int y1) {
  std::vector<std::array<int, 2>> listVoxels;
  listVoxels.push_back({x0, y0});
  const int dx= std::abs(x1 - x0);
  const int dy= std::abs(y1 - y0);
  const int xs= (x1 > x0) ? 1 : -1;
  const int ys= (y1 > y0) ? 1 : -1;

  // Driving axis is X-axis
  if (dx >= dy) {
    int p1= 2 * dy - dx;
    while (x0 != x1) {
      x0+= xs;
      if (p1 >= 0) {
        y0+= ys;
        p1-= 2 * dx;
      }
      p1+= 2 * dy;
      listVoxels.push_back({x0, y0});
    }
  }

  // Driving axis is Y-axis
  else {
    int p1= 2 * dx - dy;
    while (y0 != y1) {
      y0+= ys;
      if (p1 >= 0) {
        x0+= xs;
        p1-= 2 * dy;
      }
      p1+= 2 * dx;
      listVoxels.push_back({x0, y0});
    }
  }

  return listVoxels;
}


std::vector<std::array<int, 3>> Bresenham::Line3D(int x0, int y0, int z0, int x1, int y1, int z1) {
  std::vector<std::array<int, 3>> listVoxels;
  listVoxels.push_back({x0, y0, z0});
  const int dx= std::abs(x1 - x0);
  const int dy= std::abs(y1 - y0);
  const int dz= std::abs(z1 - z0);
  const int xs= (x1 > x0) ? 1 : -1;
  const int ys= (y1 > y0) ? 1 : -1;
  const int zs= (z1 > z0) ? 1 : -1;

  // Driving axis is X-axis
  if (dx >= dy && dx >= dz) {
    int p1= 2 * dy - dx;
    int p2= 2 * dz - dx;
    while (x0 != x1) {
      x0+= xs;
      if (p1 >= 0) {
        y0+= ys;
        p1-= 2 * dx;
      }
      if (p2 >= 0) {
        z0+= zs;
        p2-= 2 * dx;
      }
      p1+= 2 * dy;
      p2+= 2 * dz;
      listVoxels.push_back({x0, y0, z0});
    }
  }

  // Driving axis is Y-axis
  else if (dy >= dx && dy >= dz) {
    int p1= 2 * dx - dy;
    int p2= 2 * dz - dy;
    while (y0 != y1) {
      y0+= ys;
      if (p1 >= 0) {
        x0+= xs;
        p1-= 2 * dy;
      }
      if (p2 >= 0) {
        z0+= zs;
        p2-= 2 * dy;
      }
      p1+= 2 * dx;
      p2+= 2 * dz;
      listVoxels.push_back({x0, y0, z0});
    }
  }

  // Driving axis is Z-axis
  else {
    int p1= 2 * dy - dz;
    int p2= 2 * dx - dz;
    while (z0 != z1) {
      z0+= zs;
      if (p1 >= 0) {
        y0+= ys;
        p1-= 2 * dz;
      }
      if (p2 >= 0) {
        x0+= xs;
        p2-= 2 * dz;
      }
      p1+= 2 * dy;
      p2+= 2 * dx;
      listVoxels.push_back({x0, y0, z0});
    }
  }

  return listVoxels;
}
