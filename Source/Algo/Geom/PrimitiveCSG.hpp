#pragma once

// Standard lib
#include <array>
#include <vector>


// Generate and combine primitives via CSG operations (distance is negative inside and positive outside by convention)
class PrimitiveCSG
{
  public:
  enum class BooleanMode
  {
    Union,
    Intersection,
    Difference,
    DifferenceFlipOrder
  };


  static void Sphere(
      const int nX,
      const int nY,
      const int nZ,
      std::array<double, 3> const& iCenter,
      double const& iRadius,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  static void Cylinder(
      const int nX,
      const int nY,
      const int nZ,
      std::array<double, 3> const& iPointA,
      std::array<double, 3> const& iPointB,
      double const& iRadiusA,
      double const& iRadiusB,
      bool const& iUseRoundedEnds,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  static void ConeRound(
      const int nX,
      const int nY,
      const int nZ,
      std::array<double, 3> const& iPointA,
      std::array<double, 3> const& iPointB,
      double const& iRadiusA,
      double const& iRadiusB,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  static void AxisAlignedTorus(
      const int nX,
      const int nY,
      const int nZ,
      const int iDir,
      std::array<double, 3> const& iCenter,
      double const& iRadiusA,
      double const& iRadiusB,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  static void AxisAlignedBox(
      const int nX,
      const int nY,
      const int nZ,
      std::array<double, 3> const& iCornerMin,
      std::array<double, 3> const& iCornerMax,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  // Compute the distance field to the closed polyline
  // Then replicate along the axis-aligend extrusion direction until reaching the chosen limits
  static void AxisAlignedExtrudedSketch(
      const int nX,
      const int nY,
      const int nZ,
      int const iExtrusionDir,
      double const iExtruLimitMin,
      double const iExtruLimitMax,
      std::vector<std::array<double, 3>> const& iVertices,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  // Translate and rotate the field into the sketch reference frame defined by its center and vector
  // Then compute the signed distance to that sketch extruded up to the chosen limits
  static void GeneralExtrudedSketch(
      const int nX,
      const int nY,
      const int nZ,
      std::array<double, 3> const& iSketchCenter,
      std::array<double, 3> const& iSketchVector,
      double const iExtruLimitNega,
      double const iExtruLimitPosi,
      std::vector<std::array<double, 3>> const& iVertices,
      PrimitiveCSG::BooleanMode const& iMode,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<double>& ioDistanceField);


  // Boolean operation between two distance fields
  static void BooleanOperation(
      const int nX,
      const int nY,
      const int nZ,
      PrimitiveCSG::BooleanMode const& iMode,
      std::vector<double> const& iDistanceField,
      std::vector<double>& ioDistanceField);
};
