#pragma once

// Standard lib
#include <array>
#include <string>
#include <vector>


class FileInput
{
  public:
  static void LoadScalarFieldBinaryFile(
      std::string const iFullpath,
      int const iNbX,
      int const iNbY,
      int const iNbZ,
      std::vector<std::vector<std::vector<double>>>& oField,
      bool const iVerbose);

  static void LoadVectorFieldBinaryFile(
      std::string const iFullpath,
      int const iNbX,
      int const iNbY,
      int const iNbZ,
      std::vector<std::vector<std::vector<std::array<double, 3>>>>& oField,
      bool const iVerbose);

  static void LoadBoxTXTFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      bool const iVerbose);

  static void LoadScalarFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>>& oField,
      bool const iVerbose);

  static void LoadScalarFieldAsciiVTIFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      std::vector<std::vector<std::vector<double>>>& oField,
      bool const iVerbose);

  static void LoadScalarFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      std::vector<std::vector<std::vector<double>>>& oField,
      bool const iVerbose);

  static void LoadVectorFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      std::vector<std::vector<std::vector<std::array<double, 3>>>>& oField,
      bool const iVerbose);

  static void LoadScalarFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<int>>>& oField,
      bool const iVerbose);

  static void LoadImageBMPFile(
      std::string const iFullpath,
      std::vector<std::vector<std::array<float, 4>>>& oImageRGBA,
      bool const iVerbose);

  static void LoadVectorFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>>& oField,
      bool const iVerbose);

  static void LoadVectorFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 3>>>>& oField,
      bool const iVerbose);

  static void LoadTensorFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 9>>>>& oField,
      bool const iVerbose);

  static void LoadMeshOBJFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>>& oPoints,
      std::vector<std::array<double, 3>>& oColors,
      std::vector<std::array<int, 3>>& oTriangles,
      bool const iVerbose);
};
