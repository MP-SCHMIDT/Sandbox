
#include "FileOutput.hpp"

// Standard lib
#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

// Algo headers
#include "Geom/BoxGrid.hpp"
#include "Type/Field.hpp"
#include "Type/Vec.hpp"


bool FileOutput::SaveBoxTXTFile(
    std::string const iFullpath,
    int const iNbX,
    int const iNbY,
    int const iNbZ,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    bool const iVerbose) {
  if (iVerbose) printf("Saving TXT box file [%s] ", iFullpath.c_str());

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  double voxSizeX, voxSizeY, voxSizeZ, voxSizeDiag;
  BoxGrid::GetVoxelSizes(iNbX, iNbY, iNbZ, iBBoxMin, iBBoxMax, true, voxSizeX, voxSizeY, voxSizeZ, voxSizeDiag);

  fprintf(outputFile, "%f %f %f\n", iBBoxMin[0], iBBoxMin[1], iBBoxMin[2]);
  fprintf(outputFile, "%f %f %f\n", iBBoxMax[0], iBBoxMax[1], iBBoxMax[2]);
  fprintf(outputFile, "%d %d %d\n", iNbX, iNbY, iNbZ);
  fprintf(outputFile, "%f %f %f\n", voxSizeX, voxSizeY, voxSizeZ);

  fclose(outputFile);

  if (iVerbose) printf("File saved: %f x %f x %f -> %f x %f x %f \n", iBBoxMin[0], iBBoxMin[1], iBBoxMin[2], iBBoxMax[0], iBBoxMax[1], iBBoxMax[2]);
  return true;
}


bool FileOutput::SaveScalarFieldFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<double>>> const& iDoubleField,
    std::vector<std::vector<std::vector<int>>> const& iIntField,
    std::vector<std::vector<std::vector<bool>>> const& iBoolField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving TXT scalar field file [%s] ", iFullpath.c_str());

  int type= 0, nbX= 0, nbY= 0, nbZ= 0;
  if (!iDoubleField.empty() && !iDoubleField[0].empty() && !iDoubleField[0][0].empty()) {
    type= 1;
    nbX= (int)iDoubleField.size();
    nbY= (int)iDoubleField[0].size();
    nbZ= (int)iDoubleField[0][0].size();
  }
  else if (!iIntField.empty() && !iIntField[0].empty() && !iIntField[0][0].empty()) {
    type= 2;
    nbX= (int)iIntField.size();
    nbY= (int)iIntField[0].size();
    nbZ= (int)iIntField[0][0].size();
  }
  else if (!iBoolField.empty() && !iBoolField[0].empty() && !iBoolField[0][0].empty()) {
    type= 3;
    nbX= (int)iBoolField.size();
    nbY= (int)iBoolField[0].size();
    nbZ= (int)iBoolField[0][0].size();
  }

  if (type == 0) {
    printf("[ERROR] Invalid field dimensions\n");
    return false;
  }

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  fprintf(outputFile, "%d %d %d\n", nbX, nbY, nbZ);
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (type == 1) {
          if (iDoubleField[x][y][z] == 0.0)
            fprintf(outputFile, "0\n");
          else if (iDoubleField[x][y][z] == 1.0)
            fprintf(outputFile, "1\n");
          else
            fprintf(outputFile, "%.3e\n", iDoubleField[x][y][z]);
        }
        if (type == 2) fprintf(outputFile, "%d\n", iIntField[x][y][z]);
        if (type == 3) fprintf(outputFile, "%d\n", (iBoolField[x][y][z]) ? 1 : 0);
      }
    }
  }
  fclose(outputFile);

  if (iVerbose) printf("File saved: %d x %d x %d voxels\n", nbX, nbY, nbZ);
  return true;
}


bool FileOutput::SaveScalarListTXTFile(
    std::string const iFullpath,
    std::vector<double> const& iVector,
    bool const iVerbose) {
  if (iVerbose) printf("Saving vector of scalars TXT file [%s] ", iFullpath.c_str());

  if (iVector.empty()) {
    printf("[ERROR] Invalid vector dimensions\n");
    return false;
  }

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  for (int k= 0; k < int(iVector.size()); k++) {
    fprintf(outputFile, "%f\n", iVector[k]);
  }

  fclose(outputFile);

  if (iVerbose) printf("File saved: %d values\n", int(iVector.size()));
  return true;
}


bool FileOutput::SaveScalarListBinaryFile(
    std::string const iFullpath,
    std::vector<float> const& iField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving RAW scalar field binary file [%s] ", iFullpath.c_str());

  if (iField.empty()) {
    printf("[ERROR] Invalid field dimensions\n");
    return false;
  }

  std::ofstream outputFile;
  outputFile.open(iFullpath, std::ios::binary);
  if (!outputFile.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  for (int k= 0; k < (int)iField.size(); k++) {
    outputFile.write((char*)&iField[k], sizeof(float));
  }

  outputFile.close();

  if (iVerbose) printf("File saved: %d values\n", (int)iField.size());
  return true;
}


bool FileOutput::SaveSparseMatrixMarketFile(
    std::string const iFullpath,
    int const iNbRows,
    int const iNbCols,
    std::vector<int> const& iRow,
    std::vector<int> const& iCol,
    std::vector<double> const& iVal,
    bool const iVerbose) {
  if (iVerbose) printf("Saving Sparse MatrixMarket file [%s] ", iFullpath.c_str());

  if (iRow.size() != iCol.size() || iRow.size() != iVal.size()) {
    printf("[ERROR] Invalid matrix dimensions\n");
    return false;
  }

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  fprintf(outputFile, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(outputFile, "%d %d %d\n", iNbRows, iNbCols, (int)iVal.size());

  for (int k= 0; k < (int)iVal.size(); k++) {
    fprintf(outputFile, "%d %d %e\n", iRow[k] + 1, iCol[k] + 1, iVal[k]);
  }

  fclose(outputFile);

  if (iVerbose) printf("File saved: Matrix %d x %d with %d non zeros\n", iNbRows, iNbCols, (int)iVal.size());
  return true;
}


bool FileOutput::SaveVectorFieldFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iDoubleField,
    std::vector<std::vector<std::vector<std::array<int, 3>>>> const& iIntField,
    std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iBoolField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving TXT vector field file [%s] ", iFullpath.c_str());

  int type= 0, nbX= 0, nbY= 0, nbZ= 0;
  if (!iDoubleField.empty() && !iDoubleField[0].empty() && !iDoubleField[0][0].empty()) {
    type= 1;
    nbX= (int)iDoubleField.size();
    nbY= (int)iDoubleField[0].size();
    nbZ= (int)iDoubleField[0][0].size();
  }
  else if (!iIntField.empty() && !iIntField[0].empty() && !iIntField[0][0].empty()) {
    type= 2;
    nbX= (int)iIntField.size();
    nbY= (int)iIntField[0].size();
    nbZ= (int)iIntField[0][0].size();
  }
  else if (!iBoolField.empty() && !iBoolField[0].empty() && !iBoolField[0][0].empty()) {
    type= 3;
    nbX= (int)iBoolField.size();
    nbY= (int)iBoolField[0].size();
    nbZ= (int)iBoolField[0][0].size();
  }

  if (type == 0) {
    printf("[ERROR] Invalid field dimensions\n");
    return false;
  }

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  fprintf(outputFile, "%d %d %d\n", nbX, nbY, nbZ);
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (type == 1) {
          if (iDoubleField[x][y][z][0] != 0.0 || iDoubleField[x][y][z][1] != 0.0 || iDoubleField[x][y][z][2] != 0.0)
            fprintf(outputFile, "%.3e %.3e %.3e\n", iDoubleField[x][y][z][0], iDoubleField[x][y][z][1], iDoubleField[x][y][z][2]);
          else
            fprintf(outputFile, "0 0 0\n");
        }
        if (type == 2) fprintf(outputFile, "%d %d %d\n", iIntField[x][y][z][0], iIntField[x][y][z][1], iIntField[x][y][z][2]);
        if (type == 3) fprintf(outputFile, "%d %d %d\n", (iBoolField[x][y][z][0]) ? 1 : 0, (iBoolField[x][y][z][1]) ? 1 : 0, (iBoolField[x][y][z][2]) ? 1 : 0);
      }
    }
  }
  fclose(outputFile);

  if (iVerbose) printf("File saved: %d x %d x %d voxels\n", nbX, nbY, nbZ);
  return true;
}


bool FileOutput::SaveTensorFieldFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<std::array<double, 9>>>> const& iDoubleField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving TXT tensor field file [%s] ", iFullpath.c_str());

  int nbX, nbY, nbZ;
  Field::GetDim(iDoubleField, nbX, nbY, nbZ);

  if (nbX <= 0 || nbY <= 0 || nbZ <= 0) {
    printf("[ERROR] Invalid field dimensions\n");
    return false;
  }

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  fprintf(outputFile, "%d %d %d\n", nbX, nbY, nbZ);
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        for (int k= 0; k < 9; k++) {
          if (iDoubleField[x][y][z][k] == 0.0)
            fprintf(outputFile, "0");
          else if (iDoubleField[x][y][z][k] == 1.0)
            fprintf(outputFile, "1");
          else
            fprintf(outputFile, "%.3e", iDoubleField[x][y][z][k]);
          if (k < 8)
            fprintf(outputFile, " ");
        }
        fprintf(outputFile, "\n");
      }
    }
  }
  fclose(outputFile);

  if (iVerbose) printf("File saved: %d x %d x %d voxels\n", nbX, nbY, nbZ);
  return true;
}


bool FileOutput::SaveScalarFieldRawVTIFile(
    std::string const iFullpath,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>> const& iField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving raw scalar VTI file [%s] ", iFullpath.c_str());

  // Create the file
  std::ofstream outputFile;
  outputFile.open(iFullpath, std::ios::binary);
  if (!outputFile.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  // Get the voxels dimensions
  int nbX, nbY, nbZ;
  Field::GetDim(iField, nbX, nbY, nbZ);
  if (nbX <= 0 || nbY <= 0 || nbZ <= 0) return false;
  double voxSizeX, voxSizeY, voxSizeZ, voxDiag, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, voxSizeX, voxSizeY, voxSizeZ, voxDiag);
  BoxGrid::GetVoxelStart(iBBoxMin, voxSizeX, voxSizeY, voxSizeZ, true, startX, startY, startZ);

  // Add the header
  outputFile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
  // outputFile << "  <ImageData WholeExtent=\"0 " + std::to_string(nbX-1)    + " 0 " + std::to_string(nbY-1)    + " 0 " + std::to_string(nbZ-1)
  //                              + "\" Origin=\"" + std::to_string(startX)   + " "   + std::to_string(startY)   + " "   + std::to_string(startZ)
  //                             + "\" Spacing=\"" + std::to_string(voxSizeX) + " "   + std::to_string(voxSizeY) + " "   + std::to_string(voxSizeZ) + "\">" << std::endl;
  outputFile << "  <ImageData WholeExtent=\"0 " + std::to_string(nbX - 1) + " 0 " + std::to_string(nbY - 1) + " 0 " + std::to_string(nbZ - 1) + "\"";
  outputFile << +" Origin=\"" + std::to_string(startX) + " " + std::to_string(startY) + " " + std::to_string(startZ) + "\"";
  outputFile << +" Spacing=\"" + std::to_string(voxSizeX) + " " + std::to_string(voxSizeY) + " " + std::to_string(voxSizeZ) + "\">" << std::endl;
  outputFile << "    <Piece Extent=\"0 " + std::to_string(nbX - 1) + " 0 " + std::to_string(nbY - 1) + " 0 " + std::to_string(nbZ - 1) + "\">" << std::endl;
  outputFile << "      <PointData>" << std::endl;
  outputFile << "        <DataArray type=\"Float32\" Name=\"scalarField\" format=\"appended\" RangeMin=\"0\" RangeMax=\"1\" offset=\"0\">" << std::endl;
  outputFile << "        </DataArray>" << std::endl;
  outputFile << "      </PointData>" << std::endl;
  outputFile << "      <CellData>" << std::endl;
  outputFile << "      </CellData>" << std::endl;
  outputFile << "    </Piece>" << std::endl;
  outputFile << "  </ImageData>" << std::endl;
  outputFile << "  <AppendedData encoding=\"raw\">" << std::endl;

  // Add the underscore delimiter followed by the number of data bytes using header_type="UInt64"
  outputFile << "   _";
  uint64_t const nbBytes= sizeof(float) * nbX * nbY * nbZ;
  outputFile.write((char*)&nbBytes, sizeof(uint64_t));

  // Add the raw data bytes in z-major x-minor order
  for (int z= 0; z < nbZ; z++) {
    for (int y= 0; y < nbY; y++) {
      for (int x= 0; x < nbX; x++) {
        float val= float(iField[x][y][z]);
        outputFile.write((char*)&val, sizeof(float));
      }
    }
  }

  // Add the closing tags
  outputFile << std::endl;
  outputFile << "  </AppendedData>" << std::endl;
  outputFile << "</VTKFile>" << std::endl;

  // Flush and exit
  outputFile.flush();
  outputFile.close();

  if (iVerbose) printf("File saved: %d x %d x %d voxels\n", nbX, nbY, nbZ);
  return true;
}


bool FileOutput::SaveVectorFieldRawVTIFile(
    std::string const iFullpath,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving raw vector VTI file [%s] ", iFullpath.c_str());

  // Create the file
  std::ofstream outputFile;
  outputFile.open(iFullpath, std::ios::binary);
  if (!outputFile.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  // Get the voxels dimensions
  int nbX, nbY, nbZ;
  Field::GetDim(iField, nbX, nbY, nbZ);
  double voxSizeX, voxSizeY, voxSizeZ, voxDiag;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, voxSizeX, voxSizeY, voxSizeZ, voxDiag);

  // Add the header
  outputFile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl;
  outputFile << "  <ImageData WholeExtent=\"0 " + std::to_string(nbX - 1) + " 0 " + std::to_string(nbY - 1) + " 0 " + std::to_string(nbZ - 1) + "\" Origin=\"" + std::to_string(iBBoxMin[0] + 0.5 * voxSizeX) + " " + std::to_string(iBBoxMin[1] + 0.5 * voxSizeY) + " " + std::to_string(iBBoxMin[2] + 0.5 * voxSizeZ) + "\" Spacing=\"" + std::to_string(voxSizeX) + " " + std::to_string(voxSizeY) + " " + std::to_string(voxSizeZ) + "\">" << std::endl;
  outputFile << "    <Piece Extent=\"0 " + std::to_string(nbX - 1) + " 0 " + std::to_string(nbY - 1) + " 0 " + std::to_string(nbZ - 1) + "\">" << std::endl;
  outputFile << "      <PointData Scalars=\"ImageFile\">" << std::endl;
  outputFile << "        <DataArray type=\"Float32\" Name=\"vectorField\" NumberOfComponents=\"3\" format=\"appended\" RangeMin=\"0\" RangeMax=\"1\" offset=\"0\">" << std::endl;
  outputFile << "          <InformationKey name=\"L2_NORM_RANGE\" location=\"vtkDataArray\" length=\"2\">" << std::endl;
  outputFile << "            <Value index=\"0\">" << std::endl;
  outputFile << "              0" << std::endl;
  outputFile << "            </Value>" << std::endl;
  outputFile << "            <Value index=\"1\">" << std::endl;
  outputFile << "              1" << std::endl;
  outputFile << "            </Value>" << std::endl;
  outputFile << "          </InformationKey>" << std::endl;
  outputFile << "        </DataArray>" << std::endl;
  outputFile << "      </PointData>" << std::endl;
  outputFile << "      <CellData>" << std::endl;
  outputFile << "      </CellData>" << std::endl;
  outputFile << "    </Piece>" << std::endl;
  outputFile << "  </ImageData>" << std::endl;
  outputFile << "  <AppendedData encoding=\"raw\">" << std::endl;

  // Add the underscore delimiter followed by the number of data bytes using header_type="UInt64"
  outputFile << "   _";
  uint64_t const nbBytes= sizeof(float) * nbX * nbY * nbZ * 3;
  outputFile.write((char*)&nbBytes, sizeof(uint64_t));

  // Add the raw data bytes in z-major x-minor order
  for (int z= 0; z < nbZ; z++) {
    for (int y= 0; y < nbY; y++) {
      for (int x= 0; x < nbX; x++) {
        float val0= float(iField[x][y][z][0]);
        float val1= float(iField[x][y][z][1]);
        float val2= float(iField[x][y][z][2]);
        outputFile.write((char*)&val0, sizeof(float));
        outputFile.write((char*)&val1, sizeof(float));
        outputFile.write((char*)&val2, sizeof(float));
      }
    }
  }

  // Add the closing tags
  outputFile << std::endl;
  outputFile << "  </AppendedData>" << std::endl;
  outputFile << "</VTKFile>" << std::endl;

  // Flush and exit
  outputFile.flush();
  outputFile.close();

  if (iVerbose) printf("File saved: %d x %d x %d voxels\n", nbX, nbY, nbZ);
  return true;
}


bool FileOutput::SaveMeshOBJFile(
    std::string const iFullpath,
    std::vector<std::array<double, 3>> const& iVertices,
    std::vector<std::array<double, 3>> const& iVerticesColors,
    std::vector<std::array<int, 3>> const& iTriangles,
    std::vector<std::array<int, 4>> const& iQuads,
    bool const iVerbose) {
  if (iVerbose) printf("Saving OBJ mesh file [%s] ", iFullpath.c_str());

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  for (unsigned int k= 0; k < iVertices.size(); k++) {
    if (k < iVerticesColors.size())
      fprintf(outputFile, "v %lf %lf %lf %lf %lf %lf\n", iVertices[k][0], iVertices[k][1], iVertices[k][2], iVerticesColors[k][0], iVerticesColors[k][1], iVerticesColors[k][2]);
    else
      fprintf(outputFile, "v %lf %lf %lf\n", iVertices[k][0], iVertices[k][1], iVertices[k][2]);
  }
  for (unsigned int k= 0; k < iTriangles.size(); k++) {
    fprintf(outputFile, "f %d %d %d\n", iTriangles[k][0] + 1, iTriangles[k][1] + 1, iTriangles[k][2] + 1);
  }
  for (unsigned int k= 0; k < iQuads.size(); k++) {
    fprintf(outputFile, "f %d %d %d %d\n", iQuads[k][0] + 1, iQuads[k][1] + 1, iQuads[k][2] + 1, iQuads[k][3] + 1);
  }

  fclose(outputFile);

  if (iVerbose) printf("File saved: %zd vertices, %zd triangles, %zd quads\n", iVertices.size(), iTriangles.size(), iQuads.size());
  return true;
}


bool FileOutput::SaveGraphINPFile(
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
    bool const iVerbose) {
  if (iVerbose) printf("Saving INP graph file [%s] ", iFullpath.c_str());

  std::ofstream ofs;
  ofs.imbue(std::locale::classic());
  ofs.open(iFullpath);
  if (!ofs.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  // Write the material
  ofs << "*MATERIAL, NAME=MAT1" << std::endl;
  ofs << "*ELASTIC, TYPE=ISOTROPIC" << std::endl;
  ofs << iYoungModulus << ", " << iPoissonRatio << std::endl;

  // Write the nodes
  ofs << "*NODE" << std::endl;
  for (int i= 0; i < (int)iNodes.size(); i++)
    ofs << i + 1 << ", " << iNodes[i][0] << ", " << iNodes[i][1] << ", " << iNodes[i][2] << std::endl;

  int numElement= 1;

  // Write the truss elements
  // if (writeBars) {
  //  for (int i= 0; i < (int)iBars.size(); i++) {
  //    ofs << "*ELEMENT, TYPE=T3D2, ELSET=BAR_" << numElement << std::endl;
  //    ofs << numElement << ", " << iBars[i][0]+1 << ", " << iBars[i][1]+1 << std::endl;
  //    ofs << "*SOLID SECTION, ELSET=BAR_" << numElement << ", MATERIAL=MAT1" << std::endl;
  //    ofs << EIGEN_PI*iBarRadii[i]*iBarRadii[i] << std::endl;
  //    numElement++;
  //  }
  //}

  // Write the beam elements
  if (iWriteBars) {
    for (int i= 0; i < (int)iBars.size(); i++) {
      ofs << "*ELEMENT, TYPE=B31, ELSET=BAR_" << numElement << std::endl;
      ofs << numElement << ", " << iBars[i][0] + 1 << ", " << iBars[i][1] + 1 << std::endl;
      ofs << "*BEAM SECTION, SECTION=CIRC, ELSET=BAR_" << numElement << ", MATERIAL=MAT1" << std::endl;
      ofs << iBarRadii[i] << std::endl;
      Vec::Vec3<double> p0(iNodes[iBars[i][0]][0], iNodes[iBars[i][0]][1], iNodes[iBars[i][0]][2]);
      Vec::Vec3<double> p1(iNodes[iBars[i][1]][0], iNodes[iBars[i][1]][1], iNodes[iBars[i][1]][2]);
      Vec::Vec3<double> U= Vec::Vec3<double>(p1 - p0).cross(Vec::Vec3<double>(0.0, 0.0, 1.0));
      if (U.norm() < 1.e-10) U= Vec::Vec3<double>(p1 - p0).cross(Vec::Vec3<double>(0.0, 1.0, 0.0));
      U.normalize();
      ofs << U[0] << ", " << U[1] << ", " << U[2] << std::endl;
      numElement++;
    }
  }

  // Write the shell elements
  if (iWriteTris) {
    for (int i= 0; i < (int)iTris.size(); i++) {
      ofs << "*ELEMENT, TYPE=S3, ELSET=TRI_" << numElement << std::endl;
      ofs << numElement << ", " << iTris[i][0] + 1 << ", " << iTris[i][1] + 1 << ", " << iTris[i][2] + 1 << std::endl;
      ofs << "*SHELL SECTION, ELSET=TRI_" << numElement << ", MATERIAL=MAT1" << std::endl;
      ofs << iTriThicknesses[i] << std::endl;
      numElement++;
    }
  }

  // Write the volume elements
  if (iWriteTets) {
    for (int i= 0; i < (int)iTets.size(); i++) {
      ofs << "*ELEMENT, TYPE=C3D4, ELSET=TET_" << numElement << std::endl;
      ofs << numElement << ", " << iTets[i][0] + 1 << ", " << iTets[i][1] + 1 << ", " << iTets[i][2] + 1 << ", " << iTets[i][3] + 1 << std::endl;
      ofs << "*SOLID SECTION, ELSET=TET_" << numElement << ", MATERIAL=MAT1" << std::endl;
      ofs << 1 << std::endl;
      numElement++;
    }
  }

  // Write the boundary conditions
  ofs << "*BOUNDARY" << std::endl;
  for (int i= 0; i < (int)iNodes.size(); i++)
    for (int j= 0; j < 3; j++)
      if (iClamps[i][j])
        ofs << i + 1 << ", " << j + 1 << ", " << j + 1 << std::endl;

  // Write the loads
  ofs << "*STEP, NAME=STEP1, PERTURBATION" << std::endl;
  ofs << "*STATIC" << std::endl;
  ofs << "*CLOAD" << std::endl;

  for (int i= 0; i < (int)iNodes.size(); i++)
    for (int j= 0; j < 3; j++)
      if (fabs(iLoads[i][j]) > 0.0)
        ofs << i + 1 << ", " << j + 1 << ", " << iLoads[i][j] << std::endl;

  // Finish
  ofs << "*OUTPUT, FIELD, VARIABLE=PRESELECT" << std::endl;
  ofs << "*OUTPUT, HISTORY, VARIABLE=PRESELECT" << std::endl;
  ofs << "*END STEP" << std::endl;
  ofs.close();

  if (iVerbose) printf("File saved: %d nodes, %d elements\n", int(iNodes.size()), numElement - 1);
  return true;
}


bool FileOutput::SaveHexaMeshINPFile(
    std::string const iFullpath,
    std::array<double, 3> const& iBBoxMin,
    std::array<double, 3> const& iBBoxMax,
    std::vector<std::vector<std::vector<double>>> const& iDensityField,
    std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
    std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
    double const iDensityThreshold,
    bool const iVerbose) {
  if (iVerbose) printf("Saving hexa mesh INP file [%s] ", iFullpath.c_str());

  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetDim(iDensityField, nbX, nbY, nbZ);
  if (nbX <= 0 || nbY <= 0 || nbZ <= 0) return false;
  double stepX, stepY, stepZ, voxDiag, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ, voxDiag);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, false, startX, startY, startZ);


  // Create file in write/overwrite mode
  std::ofstream ofs;
  ofs.imbue(std::locale::classic());
  ofs.open(iFullpath);
  if (!ofs.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }


  // Find nodes to write based on density threshold and create the index mapping
  std::vector<std::vector<std::vector<int>>> idx(nbX + 1,
                                                 std::vector<std::vector<int>>(nbY + 1,
                                                                               std::vector<int>(nbZ + 1, -1)));
  int nodeNum= 1;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        int xMin= std::max(x - 1, 0), xMax= std::min(x, nbX - 1);
        int yMin= std::max(y - 1, 0), yMax= std::min(y, nbY - 1);
        int zMin= std::max(z - 1, 0), zMax= std::min(z, nbZ - 1);
        if (iDensityField[xMin][yMin][zMin] < iDensityThreshold && iDensityField[xMin][yMin][zMax] < iDensityThreshold && iDensityField[xMin][yMax][zMin] < iDensityThreshold && iDensityField[xMin][yMax][zMax] < iDensityThreshold && iDensityField[xMax][yMin][zMin] < iDensityThreshold && iDensityField[xMax][yMin][zMax] < iDensityThreshold && iDensityField[xMax][yMax][zMin] < iDensityThreshold && iDensityField[xMax][yMax][zMax] < iDensityThreshold)
          continue;
        idx[x][y][z]= nodeNum;
        nodeNum++;
      }
    }
  }


  //// Write a default isotropic elastoplastic steel material
  // ofs << "*Material, name=MATERIAL_DEFAULT" << std::endl;
  // ofs << "*Density" << std::endl;
  // ofs << "  7850.0," << std::endl; // Material density in kg m^3
  // ofs << "*Elastic" << std::endl;
  // ofs << "  210.0e+09, 0.3" << std::endl; // Young modulus in kg m^-1 s^-2, Poisson ratio unitless
  // ofs << "*Plastic" << std::endl;
  // ofs << "  280.0e+06, 0.000" << std::endl; // Yield stress in N m^2, Plastic strain unitless
  // ofs << "  325.0e+06, 0.025" << std::endl; // Yield stress in N m^2, Plastic strain unitless
  // ofs << "  345.0e+06, 0.060" << std::endl; // Yield stress in N m^2, Plastic strain unitless
  // ofs << "  360.0e+06, 0.100" << std::endl; // Yield stress in N m^2, Plastic strain unitless
  // ofs << "  380.0e+06, 0.200" << std::endl; // Yield stress in N m^2, Plastic strain unitless

  // Write a default isotropic elastic unitary material
  ofs << "*Material, name=MATERIAL_DEFAULT" << std::endl;
  ofs << "*Density" << std::endl;
  ofs << "  1.0," << std::endl;  // Material density in kg m^3
  ofs << "*Elastic" << std::endl;
  ofs << "  1.0, 0.3" << std::endl;  // Young modulus in kg m^-1 s^-2, Poisson ratio unitless


  // Write the nodes with 3D positions in meters
  ofs << "*Node" << std::endl;
  for (int x= 0; x < nbX + 1; x++)
    for (int y= 0; y < nbY + 1; y++)
      for (int z= 0; z < nbZ + 1; z++)
        if (idx[x][y][z] != -1)
          ofs << idx[x][y][z] << ", " << (x * stepX + startX) << ", " << (y * stepY + startY) << ", " << (z * stepZ + startZ) << std::endl;


  // Write the hexahedral elements
  int elemNum= 1;
  // ofs << "*Element, type=C3D8" << std::endl;
  ofs << "*Element, type=C3D8I" << std::endl;
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (iDensityField[x][y][z] < iDensityThreshold) continue;
        int n0= idx[x + 0][y + 0][z + 0];
        int n1= idx[x + 1][y + 0][z + 0];
        int n2= idx[x + 1][y + 1][z + 0];
        int n3= idx[x + 0][y + 1][z + 0];
        int n4= idx[x + 0][y + 0][z + 1];
        int n5= idx[x + 1][y + 0][z + 1];
        int n6= idx[x + 1][y + 1][z + 1];
        int n7= idx[x + 0][y + 1][z + 1];
        ofs << elemNum << ", " << n0 << ", " << n1 << ", " << n2 << ", " << n3 << ", " << n4 << ", " << n5 << ", " << n6 << ", " << n7 << std::endl;
        elemNum++;
      }
    }
  }


  // Write the element set with all elements
  ofs << "*Elset, elset=ELEMSET_ALLELEM, generate" << std::endl;
  ofs << 1 << ", " << elemNum - 1 << "," << std::endl;


  // Write the node set with all nodess
  ofs << "*Nset, nset=NODESET_ALLNODE, generate" << std::endl;
  ofs << 1 << ", " << nodeNum - 1 << "," << std::endl;


  //// Find the node with the highest load
  // int idxNodeMaxLoad= 1;
  // double valNodeMaxLoad= 0.0;
  // for (int x= 0; x < nbX+1; x++) {
  //   for (int y= 0; y < nbY+1; y++) {
  //     for (int z= 0; z < nbZ+1; z++) {
  //       if (idx[x][y][z] != -1) {
  //         double valNode= std::abs(iLoadField[x][y][z][0])+std::abs(iLoadField[x][y][z][1])+std::abs(iLoadField[x][y][z][2]);
  //         if (valNodeMaxLoad < valNode) {
  //           valNodeMaxLoad= valNode;
  //           idxNodeMaxLoad= idx[x][y][z];
  //         }
  //       }
  //     }
  //   }
  // }
  //
  //
  //// Write the node set with the node having the highest load
  // ofs << "*Nset, nset=NODESET_MAXLOAD" << std::endl;
  // ofs << idxNodeMaxLoad << "," << std::endl;


  // Set the material for the solid section
  ofs << "*Solid Section, elset=ELEMSET_ALLELEM, material=MATERIAL_DEFAULT" << std::endl;


  // Start writing the boundary conditions
  ofs << "*Boundary" << std::endl;

  // Find if one displacement direction is locked for all nodes
  bool isDirXLockedForAll= true;
  bool isDirYLockedForAll= true;
  bool isDirZLockedForAll= true;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] != -1) {
          if (!iClampField[x][y][z][0]) isDirXLockedForAll= false;
          if (!iClampField[x][y][z][1]) isDirYLockedForAll= false;
          if (!iClampField[x][y][z][2]) isDirZLockedForAll= false;
        }
      }
    }
  }

  // Write the global boundary conditions (node set, first locked dof, last locked dof, value of prescribed displacement)
  if (isDirXLockedForAll) ofs << "NODESET_ALLNODE, 1, 1, 0" << std::endl;
  if (isDirYLockedForAll) ofs << "NODESET_ALLNODE, 2, 2, 0" << std::endl;
  if (isDirZLockedForAll) ofs << "NODESET_ALLNODE, 3, 3, 0" << std::endl;

  // Write the individual boundary conditions (node idx, first locked dof, last locked dof, value of prescribed displacement)
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] == -1)
          continue;
        if (!isDirXLockedForAll && iClampField[x][y][z][0]) ofs << idx[x][y][z] << ", 1, 1, 0" << std::endl;
        if (!isDirYLockedForAll && iClampField[x][y][z][1]) ofs << idx[x][y][z] << ", 2, 2, 0" << std::endl;
        if (!isDirZLockedForAll && iClampField[x][y][z][2]) ofs << idx[x][y][z] << ", 3, 3, 0" << std::endl;
      }
    }
  }


  // Write the simulation step setup
  ofs << "*STEP, NAME=STEP1, PERTURBATION" << std::endl;
  ofs << "*STATIC" << std::endl;

  //// Write the simulation step setup
  // ofs << "*Step, name=STEP1, nlgeom=YES, inc=200" << std::endl;
  // ofs << "*Static, riks"                          << std::endl;
  // ofs << "1., 1., 1e-05, 1e36, 1e35, "            << std::endl;


  // Write the concentrated loads in Newtons
  ofs << "*Cload" << std::endl;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] == -1)
          continue;
        for (int k= 0; k < 3; k++)
          if (std::abs(iLoadField[x][y][z][k]) > 0.0)
            ofs << idx[x][y][z] << ", " << k + 1 << ", " << iLoadField[x][y][z][k] << std::endl;
      }
    }
  }


  //// Write the output history requests
  // ofs << "** "                                  << std::endl;
  // ofs << "** OUTPUT REQUESTS"                   << std::endl;
  // ofs << "** "                                  << std::endl;
  // ofs << "*Restart, write, frequency=0"         << std::endl;
  // ofs << "** "                                  << std::endl;
  // ofs << "** FIELD OUTPUT: F-Output-1"          << std::endl;
  // ofs << "** "                                  << std::endl;
  // ofs << "*Output, field, variable=PRESELECT"   << std::endl;
  // ofs << "** "                                  << std::endl;
  // ofs << "** HISTORY OUTPUT: H-Output-1"        << std::endl;
  // ofs << "** "                                  << std::endl;
  // ofs << "*Output, history"                     << std::endl;
  // ofs << "*Node Output, nset=NODESET_MAXLOAD"   << std::endl;
  // ofs << "CF1, CF2, CF3, U1, U2, U3"            << std::endl;


  ofs.close();

  if (iVerbose) printf("File saved: %d nodes, %d elements\n", nodeNum - 1, elemNum - 1);
  return true;
}


bool FileOutput::SaveHexaMeshElemValuesINPFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<double>>> const& iElemValueField,
    std::vector<std::vector<std::vector<int>>> const& iDesignSpaceField,
    bool const iVerbose) {
  if (iVerbose) printf("Saving hexa mesh densities INP file [%s] ", iFullpath.c_str());
  // Get field dimensions
  int nbX, nbY, nbZ;
  Field::GetDim(iElemValueField, nbX, nbY, nbZ);
  if (nbX <= 0 || nbY <= 0 || nbZ <= 0) return false;


  // Create file in write/overwrite mode
  std::ofstream ofs;
  ofs.imbue(std::locale::classic());
  ofs.open(iFullpath);
  if (!ofs.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }


  // Write the hexahedral elements
  int elemNum= 1;
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (iDesignSpaceField[x][y][z] < 0) continue;
        ofs << elemNum << " " << iElemValueField[x][y][z] << std::endl;
        elemNum++;
      }
    }
  }


  ofs.close();

  if (iVerbose) printf("File saved: %d elements\n", elemNum - 1);
  return true;
}


bool FileOutput::SaveHexaMeshWithElsetINPFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<double>>> const& iDensityField,
    std::vector<std::vector<std::vector<int>>> const& iDesignSpaceField,
    std::vector<std::vector<std::vector<std::array<bool, 3>>>> const& iClampField,
    std::vector<std::vector<std::vector<std::array<double, 3>>>> const& iLoadField,
    double const iYoungModulus,
    double const iPoissonRatio,
    double const iDensityThreshold,
    bool const iVerbose) {
  if (iVerbose) printf("Saving hexa mesh elset INP file [%s] ", iFullpath.c_str());

  std::ofstream ofs;
  ofs.imbue(std::locale::classic());
  ofs.open(iFullpath);
  if (!ofs.is_open()) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  int nbX= int(iDensityField.size());
  int nbY= int(iDensityField[0].size());
  int nbZ= int(iDensityField[0][0].size());

  // Write the material
  ofs << "*MATERIAL, NAME=DEFAULTMATERIAL" << std::endl;
  ofs << "*ELASTIC, TYPE=ISOTROPIC" << std::endl;
  ofs << iYoungModulus << ", " << iPoissonRatio << std::endl;

  ofs << "*MATERIAL, NAME=DEFAULTMATERIALFROZEN" << std::endl;
  ofs << "*ELASTIC, TYPE=ISOTROPIC" << std::endl;
  ofs << iYoungModulus << ", " << iPoissonRatio << std::endl;

  // Find nodes to write
  std::vector<std::vector<std::vector<int>>> idx(nbX + 1,
                                                 std::vector<std::vector<int>>(nbY + 1,
                                                                               std::vector<int>(nbZ + 1, -1)));
  int nodeNum= 1;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        int xMin= std::max(x - 1, 0), xMax= std::min(x, nbX - 1);
        int yMin= std::max(y - 1, 0), yMax= std::min(y, nbY - 1);
        int zMin= std::max(z - 1, 0), zMax= std::min(z, nbZ - 1);
        if (iDensityField[xMin][yMin][zMin] < iDensityThreshold && iDensityField[xMin][yMin][zMax] < iDensityThreshold && iDensityField[xMin][yMax][zMin] < iDensityThreshold && iDensityField[xMin][yMax][zMax] < iDensityThreshold && iDensityField[xMax][yMin][zMin] < iDensityThreshold && iDensityField[xMax][yMin][zMax] < iDensityThreshold && iDensityField[xMax][yMax][zMin] < iDensityThreshold && iDensityField[xMax][yMax][zMax] < iDensityThreshold)
          continue;
        idx[x][y][z]= nodeNum;
        nodeNum++;
      }
    }
  }

  // Write the nodes
  ofs << "*NODE" << std::endl;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] != -1) {
          ofs << idx[x][y][z] << ", " << x << ", " << y << ", " << z << std::endl;
        }
      }
    }
  }

  // Write the hexahedral elements
  int numElement= 1;
  ofs << "*ELEMENT, TYPE=C3D8, ELSET=HEXAMESH" << std::endl;
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (iDensityField[x][y][z] < iDensityThreshold) continue;
        if (iDesignSpaceField[x][y][z] != 0) continue;
        int n0= idx[x + 0][y + 0][z + 0];
        int n1= idx[x + 1][y + 0][z + 0];
        int n2= idx[x + 1][y + 1][z + 0];
        int n3= idx[x + 0][y + 1][z + 0];
        int n4= idx[x + 0][y + 0][z + 1];
        int n5= idx[x + 1][y + 0][z + 1];
        int n6= idx[x + 1][y + 1][z + 1];
        int n7= idx[x + 0][y + 1][z + 1];
        ofs << numElement << ", " << n0 << ", " << n1 << ", " << n2 << ", " << n3 << ", " << n4 << ", " << n5 << ", " << n6 << ", " << n7 << std::endl;
        numElement++;
      }
    }
  }
  ofs << "*ELEMENT, TYPE=C3D8, ELSET=HEXAMESHFROZEN" << std::endl;
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (iDensityField[x][y][z] < iDensityThreshold) continue;
        if (iDesignSpaceField[x][y][z] == 0) continue;
        int n0= idx[x + 0][y + 0][z + 0];
        int n1= idx[x + 1][y + 0][z + 0];
        int n2= idx[x + 1][y + 1][z + 0];
        int n3= idx[x + 0][y + 1][z + 0];
        int n4= idx[x + 0][y + 0][z + 1];
        int n5= idx[x + 1][y + 0][z + 1];
        int n6= idx[x + 1][y + 1][z + 1];
        int n7= idx[x + 0][y + 1][z + 1];
        ofs << numElement << ", " << n0 << ", " << n1 << ", " << n2 << ", " << n3 << ", " << n4 << ", " << n5 << ", " << n6 << ", " << n7 << std::endl;
        numElement++;
      }
    }
  }

  // Write the material section
  ofs << "*SOLID SECTION, ELSET=HEXAMESH, MATERIAL=DEFAULTMATERIAL" << std::endl;
  ofs << "*SOLID SECTION, ELSET=HEXAMESHFROZEN, MATERIAL=DEFAULTMATERIALFROZEN" << std::endl;

  // Write the boundary conditions
  ofs << "*BOUNDARY" << std::endl;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] == -1)
          continue;
        for (int k= 0; k < 3; k++)
          if (iClampField[x][y][z][k])
            ofs << idx[x][y][z] << ", " << k + 1 << ", " << k + 1 << std::endl;
      }
    }
  }

  // Write the loads
  ofs << "*STEP, NAME=STEP1, PERTURBATION" << std::endl;
  ofs << "*STATIC" << std::endl;
  ofs << "*CLOAD" << std::endl;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] == -1)
          continue;
        for (int k= 0; k < 3; k++)
          if (std::abs(iLoadField[x][y][z][k]) > 0.0)
            ofs << idx[x][y][z] << ", " << k + 1 << ", " << iLoadField[x][y][z][k] << std::endl;
      }
    }
  }

  ofs.close();

  if (iVerbose) printf("File saved: %d nodes, %d elements\n", nodeNum - 1, numElement - 1);
  return true;
}


bool FileOutput::SaveOrthotropicHexaMeshINPFile(
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
    bool const iVerbose) {
  if (iVerbose) printf("Saving hexa mesh ortho elem INP file [%s] ", iFullpath.c_str());
  int nbX, nbY, nbZ;
  Field::GetDim(iDensityField, nbX, nbY, nbZ);
  if (nbX == 0 || nbY == 0 || nbZ == 0) return false;
  double stepX, stepY, stepZ, voxDiag, startX, startY, startZ;
  BoxGrid::GetVoxelSizes(nbX, nbY, nbZ, iBBoxMin, iBBoxMax, true, stepX, stepY, stepZ, voxDiag);
  BoxGrid::GetVoxelStart(iBBoxMin, stepX, stepY, stepZ, false, startX, startY, startZ);

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  // Write the material
  fprintf(outputFile, "*MATERIAL, NAME=DEFAULTMATERIAL\n");
  fprintf(outputFile, "*ELASTIC, TYPE=ENGINEERING CONSTANTS\n");
  fprintf(outputFile, "%f, %f, %f, %f, %f, %f, %f, %f\n", iYoungX, iYoungY, iYoungZ, iPoissonXY, iPoissonXZ, iPoissonYZ, iShearXY, iShearXZ);
  fprintf(outputFile, "%f\n", iShearYZ);

  // Find nodes to write
  std::vector<std::vector<std::vector<int>>> idx(nbX + 1,
                                                 std::vector<std::vector<int>>(nbY + 1,
                                                                               std::vector<int>(nbZ + 1, -1)));
  int nodeNum= 1;
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        int xMin= std::max(x - 1, 0), xMax= std::min(x, nbX - 1);
        int yMin= std::max(y - 1, 0), yMax= std::min(y, nbY - 1);
        int zMin= std::max(z - 1, 0), zMax= std::min(z, nbZ - 1);
        if (iDensityField[xMin][yMin][zMin] < iDensityThreshold && iDensityField[xMin][yMin][zMax] < iDensityThreshold && iDensityField[xMin][yMax][zMin] < iDensityThreshold && iDensityField[xMin][yMax][zMax] < iDensityThreshold && iDensityField[xMax][yMin][zMin] < iDensityThreshold && iDensityField[xMax][yMin][zMax] < iDensityThreshold && iDensityField[xMax][yMax][zMin] < iDensityThreshold && iDensityField[xMax][yMax][zMax] < iDensityThreshold)
          continue;
        idx[x][y][z]= nodeNum;
        nodeNum++;
      }
    }
  }

  // Write the nodes
  fprintf(outputFile, "*NODE\n");
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] != -1) {
          fprintf(outputFile, "%d, %f, %f, %f\n", idx[x][y][z], x * stepX + startX, y * stepY + startY, z * stepZ + startZ);
        }
      }
    }
  }

  // Write the hexahedral elements
  int numElement= 1;
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (iDensityField[x][y][z] < iDensityThreshold) continue;
        int n0= idx[x + 0][y + 0][z + 0];
        int n1= idx[x + 1][y + 0][z + 0];
        int n2= idx[x + 1][y + 1][z + 0];
        int n3= idx[x + 0][y + 1][z + 0];
        int n4= idx[x + 0][y + 0][z + 1];
        int n5= idx[x + 1][y + 0][z + 1];
        int n6= idx[x + 1][y + 1][z + 1];
        int n7= idx[x + 0][y + 1][z + 1];
        // Build orthogonal basis
        Vec::Vec3<double> yAxis= Vec::Vec3<double>(iFiberOrientationField[x][y][z][0], iFiberOrientationField[x][y][z][1], iFiberOrientationField[x][y][z][2]);
        Vec::Vec3<double> zAxis= (yAxis).cross(Vec::Vec3<double>(1.0, 0.0, 1.0));
        if (zAxis.norm() < 1.e-6 * yAxis.norm())
          zAxis= (yAxis).cross(Vec::Vec3<double>(0.0, 1.0, 0.0));
        zAxis.normalize();
        Vec::Vec3<double> xAxis= (yAxis).cross(zAxis);
        xAxis.normalize();
        // Write values
        fprintf(outputFile, "*ELEMENT, TYPE=C3D8, ELSET=HEXAMESH-ELEM%d\n", numElement);
        fprintf(outputFile, "%d, %d, %d, %d, %d, %d, %d, %d, %d\n", numElement, n0, n1, n2, n3, n4, n5, n6, n7);
        fprintf(outputFile, "*Orientation, name=HEXAMESH-ORIENT%d\n", numElement);
        fprintf(outputFile, "%f, %f, %f, %f, %f, %f, %f, %f, %f\n", xAxis[0], xAxis[1], xAxis[2], yAxis[0], yAxis[1], yAxis[2], x * stepX + 0.5 * stepX + startX, y * stepY + 0.5 * stepY + startY, z * stepZ + 0.5 * stepZ + startZ);
        fprintf(outputFile, "*Solid Section, elset=HEXAMESH-ELEM%d, orientation=HEXAMESH-ORIENT%d, stack direction=1, material=DEFAULTMATERIAL\n", numElement, numElement);
        numElement++;
      }
    }
  }

  // Write the boundary conditions
  fprintf(outputFile, "*BOUNDARY\n");
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] == -1)
          continue;
        for (int k= 0; k < 3; k++)
          if (iClampField[x][y][z][k])
            fprintf(outputFile, "%d, %d, %d\n", idx[x][y][z], k + 1, k + 1);
      }
    }
  }

  // Write the loads
  fprintf(outputFile, "*STEP, NAME=STEP1, PERTURBATION\n");
  fprintf(outputFile, "*STATIC\n");
  fprintf(outputFile, "*CLOAD\n");
  for (int x= 0; x < nbX + 1; x++) {
    for (int y= 0; y < nbY + 1; y++) {
      for (int z= 0; z < nbZ + 1; z++) {
        if (idx[x][y][z] == -1)
          continue;
        for (int k= 0; k < 3; k++)
          if (std::abs(iLoadField[x][y][z][k]) > 0.0)
            fprintf(outputFile, "%d, %d, %f\n", idx[x][y][z], k + 1, iLoadField[x][y][z][k]);
      }
    }
  }

  fclose(outputFile);

  if (iVerbose) printf("File saved: %d nodes, %d elements\n", nodeNum - 1, numElement - 1);
  return true;
}


bool FileOutput::SaveGraphO3PFile(
    std::string const iFullpath,
    std::vector<std::array<double, 3>> const& iNodes,
    std::vector<std::array<int, 2>> const& iBars,
    std::vector<std::array<int, 3>> const& iTris,
    std::vector<std::array<int, 4>> const& iTets,
    std::vector<double> const& iBarRadii,
    bool const iWriteBars,
    bool const iWriteTris,
    bool const iWriteTets,
    bool const iVerbose) {
  if (iVerbose) printf("Saving O3P graph file [%s] ", iFullpath.c_str());

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  for (unsigned int k= 0; k < iNodes.size(); k++)
    fprintf(outputFile, "n %f %f %f\n", iNodes[k][0], iNodes[k][1], iNodes[k][2]);

  if (iWriteBars)
    for (unsigned int k= 0; k < iBars.size(); k++)
      fprintf(outputFile, "b %d %d %f\n", iBars[k][0], iBars[k][1], iBarRadii[k]);

  if (iWriteTris)
    for (unsigned int k= 0; k < iTris.size(); k++)
      fprintf(outputFile, "f %d %d %d\n", iTris[k][0], iTris[k][1], iTris[k][2]);

  if (iWriteTets)
    for (unsigned int k= 0; k < iTets.size(); k++)
      fprintf(outputFile, "v %d %d %d %d\n", iTets[k][0], iTets[k][1], iTets[k][2], iTets[k][3]);

  fclose(outputFile);
  if (iVerbose) printf("File saved: %zd nodes, %zd bars\n", iNodes.size(), iBars.size());
  return true;
}


bool FileOutput::SaveGraphTXTFile(
    std::string const iFullpath,
    std::vector<std::array<double, 3>> const& iNodes,
    std::vector<std::array<double, 3>> const& iLoads,
    std::vector<std::array<int, 3>> const& iClamps,
    std::vector<std::array<int, 2>> const& iBars,
    std::vector<double> const& iBarRadii,
    bool const iVerbose) {
  if (iVerbose) printf("Saving TXT graph file [%s] ", iFullpath.c_str());

  if (iNodes.size() != iLoads.size()) return false;
  if (iNodes.size() != iClamps.size()) return false;

  FILE* outputFile= nullptr;
  outputFile= fopen(iFullpath.c_str(), "w");
  if (outputFile == nullptr) {
    printf("[ERROR] Unable to create the file\n");
    return false;
  }

  fprintf(outputFile, "%d %d\n", (int)iNodes.size(), (int)iBars.size());

  for (unsigned int k= 0; k < iNodes.size(); k++)
    fprintf(outputFile, "n %f %f %f %f %f %f %d %d %d\n",
            iNodes[k][0], iNodes[k][1], iNodes[k][2],
            iLoads[k][0], iLoads[k][1], iLoads[k][2],
            iClamps[k][0], iClamps[k][1], iClamps[k][2]);

  for (unsigned int k= 0; k < iBars.size(); k++)
    fprintf(outputFile, "b %d %d %f\n", iBars[k][0], iBars[k][1], iBarRadii[k]);

  fclose(outputFile);
  if (iVerbose) printf("File saved: %zd nodes, %zd bars\n", iNodes.size(), iBars.size());
  return true;
}
