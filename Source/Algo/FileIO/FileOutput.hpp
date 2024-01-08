#pragma once

// Standard lib
#include <array>
#include <vector>
#include <string>


class FileOutput
{
  public:
  static void SaveVectorBinaryFile(
      std::string const iFullpath,
      std::vector<double> const& iField,
      bool const iVerbose);

  static void SaveScalarFieldBinaryFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& field,
      bool const iVerbose);

  static void SaveVectorFieldBinaryFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& field,
      bool const iVerbose);

  static void SaveScalarFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<double>>> const& iField,
      bool const iVerbose);

  static void SaveVectorFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iField,
      bool const iVerbose);

  static void SaveBoxTXTFile(
      std::string const iFullpath,
      int const iNbX,
      int const iNbY,
      int const iNbZ,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      bool const iVerbose);

  static void SaveScalarVectorTXTFile(
      std::string const iFullpath,
      std::vector<double> const& iVector,
      bool const iVerbose);

  static void SaveScalarFieldFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& iDoubleField,
      std::vector<std::vector<std::vector<int>>> const& iIntField,
      std::vector<std::vector<std::vector<bool>>> const& iBoolField,
      bool const iVerbose);

  static void SaveVectorFieldFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iDoubleField,
      std::vector<std::vector<std::vector<std::array<int, 3>>>> const& iIntField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iBoolField,
      bool const iVerbose);

  static void SaveVectorFieldCSVFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iVectorField,
      bool const iVerbose);

  static void SaveTensorFieldFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 9>>>> const& iDoubleField,
      bool const iVerbose);

  static void SaveSparseMatrixMarketFile(
      std::string const iFullpath,
      int const iNbRows,
      int const iNbCols,
      std::vector<int> const& iRow,
      std::vector<int> const& iCol,
      std::vector<double> const& iVal,
      bool const iVerbose);

  static void SavePointCloudFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& field,
      double const cutoff,
      bool const iVerbose);

  static void SaveCATMathPointCloudFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iVertices,
      bool const iVerbose);

  static void SaveMeshOBJFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iVertices,
      std::vector<std::array<double, 3>> const& iVerticesColors,
      std::vector<std::array<int, 3>> const& iTriangles,
      std::vector<std::array<int, 4>> const& iQuads,
      bool const iVerbose);

  static void SaveGraphINPFile(
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

  static void SaveHexaMeshINPFile(
      std::string const iFullpath,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<double>>> const& iDensityField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
      double const iDensityThreshold,
      bool const iVerbose);

  static void SaveHexaMeshINPFile_ElemValues(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& iElemValueField,
      std::vector<std::vector<std::vector<int>>> const& iDesignSpaceField,
      bool const iVerbose);

  static void SaveHexaMeshWithElsetINPFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>> const& iDensityField,
      std::vector<std::vector<std::vector<int>>> const& iDesignSpaceField,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
      std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
      double const iYoungModulus,
      double const iPoissonRatio,
      double const iDensityThreshold,
      bool const iVerbose);

  static void SaveOrthotropicHexaMeshINPFile(
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

  static void SaveGraphO3PFile(
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

  static void SaveGraphTXTFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>> const& iNodes,
      std::vector<std::array<double, 3>> const& iLoads,
      std::vector<std::array<int, 3>> const& iClamps,
      std::vector<std::array<int, 2>> const& iBars,
      std::vector<double> const& iBarRadii,
      bool const iVerbose);
};
