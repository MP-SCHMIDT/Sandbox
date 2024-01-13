#pragma once

// Standard lib
#include <array>
#include <string>
#include <vector>


class FileOutput
{
  public:
  static bool SaveBoxTXTFile(
      std::string const iFullpath,
      int const iNbX,
      int const iNbY,
      int const iNbZ,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      bool const iVerbose);

  static bool SaveScalarFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<double>>> const& iField,
      bool const iVerbose);

  static bool SaveVectorFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iField,
      bool const iVerbose);

  static bool SaveScalarListTXTFile(
      std::string const iFullpath,
      std::vector<double> const& iVector,
      bool const iVerbose);

  static bool SaveScalarListBinaryFile(
      std::string const iFullpath,
      std::vector<double> const& iField,
      bool const iVerbose);

  static bool SaveScalarFieldFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& iDoubleField,
      std::vector<std::vector<std::vector<int>>> const& iIntField,
      std::vector<std::vector<std::vector<bool>>> const& iBoolField,
      bool const iVerbose);

  static bool SaveVectorFieldFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iDoubleField,
      std::vector<std::vector<std::vector<std::array<int, 3>>>> const& iIntField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iBoolField,
      bool const iVerbose);

  static bool SaveTensorFieldFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 9>>>> const& iDoubleField,
      bool const iVerbose);

  static bool SaveSparseMatrixMarketFile(
      std::string const iFullpath,
      int const iNbRows,
      int const iNbCols,
      std::vector<int> const& iRow,
      std::vector<int> const& iCol,
      std::vector<double> const& iVal,
      bool const iVerbose);

  static bool SaveMeshOBJFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iVertices,
      std::vector<std::array<double, 3>> const& iVerticesColors,
      std::vector<std::array<int, 3>> const& iTriangles,
      std::vector<std::array<int, 4>> const& iQuads,
      bool const iVerbose);

  static bool SaveGraphINPFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iNodes,
      std::vector<std::array<int, 2>> const& iBars,
      std::vector<std::array<int, 3>> const& iTris,
      std::vector<std::array<int, 4>> const& iTets,
      std::vector<std::array<int, 3>> const& iClamps,
      std::vector<std::array<double, 3>> const& iLoads,
      std::vector<double> const& iBarRadii,
      std::vector<double> const& iTriThicknesses,
      double const iYoungModulus,
      double const iPoissonRatio,
      bool const iWriteBars,
      bool const iWriteTris,
      bool const iWriteTets,
      bool const iVerbose);

  static bool SaveHexaMeshINPFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<double>>> const& iDensityField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
      double const iDensityThreshold,
      bool const iVerbose);

  static bool SaveHexaMeshElemValuesINPFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& iElemValueField,
      std::vector<std::vector<std::vector<int>>> const& iDesignSpaceField,
      bool const iVerbose);

  static bool SaveHexaMeshWithElsetINPFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& iDensityField,
      std::vector<std::vector<std::vector<int>>> const& iDesignSpaceField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
      double const iYoungModulus,
      double const iPoissonRatio,
      double const iDensityThreshold,
      bool const iVerbose);

  static bool SaveOrthotropicHexaMeshINPFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<double>>> const& iDensityField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iFiberOrientationField,
      double const iYoungX,
      double const iYoungY,
      double const iYoungZ,
      double const iPoissonXY,
      double const iPoissonXZ,
      double const iPoissonYZ,
      double const iShearXY,
      double const iShearXZ,
      double const iShearYZ,
      double const iDensityThreshold,
      bool const iVerbose);

  static bool SaveGraphO3PFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iNodes,
      std::vector<std::array<int, 2>> const& iBars,
      std::vector<std::array<int, 3>> const& iTris,
      std::vector<std::array<int, 4>> const& iTets,
      std::vector<double> const& iBarRadii,
      bool const iWriteBars,
      bool const iWriteTris,
      bool const iWriteTets,
      bool const iVerbose);

  static bool SaveGraphTXTFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iNodes,
      std::vector<std::array<double, 3>> const& iLoads,
      std::vector<std::array<int, 3>> const& iClamps,
      std::vector<std::array<int, 2>> const& iBars,
      std::vector<double> const& iBarRadii,
      bool const iVerbose);
};
