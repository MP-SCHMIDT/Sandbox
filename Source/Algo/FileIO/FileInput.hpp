#pragma once

// Standard lib
#include <array>
#include <string>
#include <vector>


class FileInput
{
  public:
  static bool LoadBoxTXTFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      bool const iVerbose);

  static bool LoadScalarListTXTFile(
      std::string const iFullpath,
      std::vector<double>& oVector,
      bool const iVerbose);

  static bool LoadScalarListBinaryFile(
      std::string const iFullpath,
      std::vector<float>& ioVector,
      bool const iVerbose);

  static bool LoadScalarFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<int>>>& oField,
      bool const iVerbose);

  static bool LoadScalarFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<double>>>& oField,
      bool const iVerbose);

  static bool LoadVectorFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<bool, 3>>>>& oField,
      bool const iVerbose);

  static bool LoadVectorFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 3>>>>& oField,
      bool const iVerbose);

  static bool LoadTensorFieldTXTFile(
      std::string const iFullpath,
      std::vector<std::vector<std::vector<std::array<double, 9>>>>& oField,
      bool const iVerbose);

  static bool LoadScalarFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      std::vector<std::vector<std::vector<double>>>& oField,
      bool const iVerbose);

  static bool LoadVectorFieldRawVTIFile(
      std::string const iFullpath,
      std::array<double, 3>& oBBoxMin,
      std::array<double, 3>& oBBoxMax,
      std::vector<std::vector<std::vector<std::array<double, 3>>>>& oField,
      bool const iVerbose);

  static bool LoadImageBMPFile(
      std::string const iFullpath,
      std::vector<std::vector<std::array<float, 4>>>& oImageRGBA,
      bool const iVerbose);

  static bool LoadMeshOBJFile(
      std::string const iFullpath,
      std::vector<std::array<double, 3>>& oPoints,
      std::vector<std::array<double, 3>>& oColors,
      std::vector<std::array<int, 3>>& oTriangles,
      bool const iVerbose);
};
