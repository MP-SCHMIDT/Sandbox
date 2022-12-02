#include "SrtFileInput.hpp"

// Standard lib
#include <fstream>
#include <iostream>

// Project
#include "SrtGenericPNGParser.hpp"


void SrtFileInput::LoadScalarFieldBinaryFile(
    std::string const iFullpath,
    int const iNbX,
    int const iNbY,
    int const iNbZ,
    std::vector<std::vector<std::vector<double>>>& oField,
    bool const iVerbose) {
  if (iNbX < 1 || iNbY < 1 || iNbZ < 1) {
    printf("[ERROR] Invalid field dimensions\n\n");
    throw 0;
  }

  if (iVerbose)
    printf("Loading RAW scalar field file [%s]\n", iFullpath.c_str());

  std::ifstream inputFile;
  inputFile.open(iFullpath, std::ios::binary);
  if (!inputFile.is_open()) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  oField=
      std::vector<std::vector<std::vector<double>>>(iNbX,
                                                    std::vector<std::vector<double>>(iNbY,
                                                                                     std::vector<double>(iNbZ)));

  for (int z= 0; z < iNbZ; z++) {
    for (int y= 0; y < iNbY; y++) {
      for (int x= 0; x < iNbX; x++) {
        inputFile.read((char*)&oField[x][y][z], sizeof(double));
      }
    }
  }

  inputFile.close();

  if (iVerbose)
    printf("File loaded: %d x %d x %d voxels\n", (int)oField.size(), (int)oField[0].size(), (int)oField[0][0].size());
}


void SrtFileInput::LoadVectorFieldBinaryFile(
    std::string const iFullpath,
    int const iNbX,
    int const iNbY,
    int const iNbZ,
    std::vector<std::vector<std::vector<Eigen::Vector3d>>>& oField,
    bool const iVerbose) {
  if (iNbX < 1 || iNbY < 1 || iNbZ < 1) {
    printf("[ERROR] Invalid field dimensions\n\n");
    throw 0;
  }

  if (iVerbose)
    printf("Loading RAW scalar field file [%s]\n", iFullpath.c_str());

  std::ifstream inputFile;
  inputFile.open(iFullpath, std::ios::binary);
  if (!inputFile.is_open()) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  oField=
      std::vector<std::vector<std::vector<Eigen::Vector3d>>>(iNbX,
                                                             std::vector<std::vector<Eigen::Vector3d>>(iNbY,
                                                                                                       std::vector<Eigen::Vector3d>(iNbZ)));

  for (int z= 0; z < iNbZ; z++) {
    for (int y= 0; y < iNbY; y++) {
      for (int x= 0; x < iNbX; x++) {
        inputFile.read((char*)&oField[x][y][z][0], sizeof(double));
        inputFile.read((char*)&oField[x][y][z][1], sizeof(double));
        inputFile.read((char*)&oField[x][y][z][2], sizeof(double));
      }
    }
  }

  inputFile.close();

  if (iVerbose)
    printf("File loaded: %d x %d x %d voxels\n", (int)oField.size(), (int)oField[0].size(), (int)oField[0][0].size());
}


void SrtFileInput::LoadBoxTXTFile(
    std::string const iFullpath,
    std::array<double, 3>& oBBoxMin,
    std::array<double, 3>& oBBoxMax,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT box file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  char buffer[1000];
  double xMin= 0.0, yMin= 0.0, zMin= 0.0;
  double xMax= 0.0, yMax= 0.0, zMax= 0.0;

  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%lf %lf %lf", &xMin, &yMin, &zMin);
  }
  else {
    printf("[ERROR] Invalid box coordinates\n");
    throw 0;
  }

  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%lf %lf %lf", &xMax, &yMax, &zMax);
  }
  else {
    printf("[ERROR] Invalid box coordinates\n");
    throw 0;
  }

  if (xMin >= xMax || yMin >= yMax || zMin >= zMax) {
    printf("[ERROR] Invalid box coordinates\n");
    throw 0;
  }

  oBBoxMin= {xMin, yMin, zMin};
  oBBoxMax= {xMax, yMax, zMax};

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %f x %f x %f -> %f x %f x %f \n", xMin, yMin, zMin, xMax, yMax, zMax);
}


void SrtFileInput::LoadScalarFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<double>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT scalar field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int X= 0, Y= 0, Z= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &X, &Y, &Z);
  }
  if (X <= 0 || Y <= 0 || Z <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<double>>>(X,
                                                        std::vector<std::vector<double>>(Y,
                                                                                         std::vector<double>(Z, 0.0)));

  // Load the field values
  for (int x= 0; x < X; x++) {
    for (int y= 0; y < Y; y++) {
      for (int z= 0; z < Z; z++) {
        double density= 0.0;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%lf", &density);
        oField[x][y][z]= density;
        // oField[x][y][z]= std::max(std::min(density, 1.0), 0.0);
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", X, Y, Z);
}


void SrtFileInput::LoadScalarFieldRawVTIFile(
    std::string const iFullpath,
    std::array<double, 3>& oBBoxMin,
    std::array<double, 3>& oBBoxMax,
    std::vector<std::vector<std::vector<double>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Reading file [%s] ", iFullpath.c_str());

  // Open the file
  std::ifstream inputFile;
  inputFile.open(iFullpath, std::ios::binary);
  if (!inputFile.is_open()) {
    if (iVerbose)
      printf("[ERROR] Unable to open the file\n");
    throw 0;
  }

  // Read the header to get dimensions and stop at the beginning of the raw data
  std::string line;
  int nbX= 0, nbY= 0, nbZ= 0;
  int readState= 0;
  while (std::getline(inputFile, line)) {
    if (readState == 0 && line.find("<ImageData") != std::string::npos) {
      int begX, begY, begZ;
      int endX, endY, endZ;
      double oriX, oriY, oriZ;
      double spaX, spaY, spaZ;
      if (std::sscanf(line.c_str(), "  <ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%lf %lf %lf\" Spacing=\"%lf %lf %lf\">",
                      &begX, &endX, &begY, &endY, &begZ, &endZ, &oriX, &oriY, &oriZ, &spaX, &spaY, &spaZ) == 12) {
        readState= 1;
        nbX= endX - begX + 1;
        nbY= endY - begY + 1;
        nbZ= endZ - begZ + 1;
        oField= std::vector<std::vector<std::vector<double>>>(nbX, std::vector<std::vector<double>>(nbY, std::vector<double>(nbZ, 0.0)));
        oBBoxMin= {oriX - 0.5 * spaX, oriY - 0.5 * spaY, oriZ - 0.5 * spaZ};
        oBBoxMax= {oriX - 0.5 * spaX + nbX * spaX, oriY - 0.5 * spaY + nbY * spaY, oriZ - 0.5 * spaZ + nbZ * spaZ};
      }
    }
    if (readState == 1 && line.find("<AppendedData encoding=\"raw\">") != std::string::npos) {
      break;
    }
  }

  // Ignore until the underscore delimiter
  while (true) {
    char c= '0';
    inputFile.read((char*)&c, sizeof(char));
    if (c == '_')
      break;
  }

  // Read the number of bytes assuming the VTI file uses header_type="UInt64"
  uint64_t nbBytes;
  inputFile.read((char*)&nbBytes, sizeof(uint64_t));

  // Read the raw data
  for (int z= 0; z < nbZ; z++) {
    for (int y= 0; y < nbY; y++) {
      for (int x= 0; x < nbX; x++) {
        float val= NAN;
        inputFile.read((char*)&val, sizeof(float));
        oField[x][y][z]= double(val);
      }
    }
  }

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", nbX, nbY, nbZ);

  inputFile.close();
}


void SrtFileInput::LoadVectorFieldRawVTIFile(
    std::string const iFullpath,
    std::array<double, 3>& oBBoxMin,
    std::array<double, 3>& oBBoxMax,
    std::vector<std::vector<std::vector<Eigen::Vector3d>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Reading file [%s] ", iFullpath.c_str());

  // Open the file
  std::ifstream inputFile;
  inputFile.open(iFullpath, std::ios::binary);
  if (!inputFile.is_open()) {
    if (iVerbose)
      printf("[ERROR] Unable to open the file\n");
    throw 0;
  }

  // Read the header to get dimensions and stop at the beginning of the raw data
  std::string line;
  int nbX= 0, nbY= 0, nbZ= 0;
  int readState= 0;
  while (std::getline(inputFile, line)) {
    if (readState == 0 && line.find("<ImageData") != std::string::npos) {
      int begX, begY, begZ;
      int endX, endY, endZ;
      double oriX, oriY, oriZ;
      double spaX, spaY, spaZ;
      if (std::sscanf(line.c_str(), "  <ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%lf %lf %lf\" Spacing=\"%lf %lf %lf\">",
                      &begX, &endX, &begY, &endY, &begZ, &endZ, &oriX, &oriY, &oriZ, &spaX, &spaY, &spaZ) == 12) {
        readState= 1;
        nbX= endX - begX + 1;
        nbY= endY - begY + 1;
        nbZ= endZ - begZ + 1;
        oField= std::vector<std::vector<std::vector<Eigen::Vector3d>>>(nbX,
                                                                       std::vector<std::vector<Eigen::Vector3d>>(nbY,
                                                                                                                 std::vector<Eigen::Vector3d>(nbZ,
                                                                                                                                              Eigen::Vector3d(0.0, 0.0, 0.0))));
        oBBoxMin= {oriX - 0.5 * spaX, oriY - 0.5 * spaY, oriZ - 0.5 * spaZ};
        oBBoxMax= {oriX - 0.5 * spaX + nbX * spaX, oriY - 0.5 * spaY + nbY * spaY, oriZ - 0.5 * spaZ + nbZ * spaZ};
      }
    }
    if (readState == 1 && line.find("<AppendedData encoding=\"raw\">") != std::string::npos) {
      break;
    }
  }

  // Ignore until the underscore delimiter
  while (true) {
    char c= '0';
    inputFile.read((char*)&c, sizeof(char));
    if (c == '_')
      break;
  }

  // Read the number of bytes assuming the VTI file uses header_type="UInt64"
  uint64_t nbBytes;
  inputFile.read((char*)&nbBytes, sizeof(uint64_t));

  // Read the raw data
  for (int z= 0; z < nbZ; z++) {
    for (int y= 0; y < nbY; y++) {
      for (int x= 0; x < nbX; x++) {
        for (int k= 0; k < 3; k++) {
          float val= NAN;
          inputFile.read((char*)&val, sizeof(float));
          oField[x][y][z][k]= double(val);
        }
      }
    }
  }

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", nbX, nbY, nbZ);

  inputFile.close();
}


void SrtFileInput::LoadScalarFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<int>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT scalar field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int X= 0, Y= 0, Z= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &X, &Y, &Z);
  }
  if (X <= 0 || Y <= 0 || Z <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<int>>>(X,
                                                     std::vector<std::vector<int>>(Y,
                                                                                   std::vector<int>(Z, 0)));

  // Load the field values
  for (int x= 0; x < X; x++) {
    for (int y= 0; y < Y; y++) {
      for (int z= 0; z < Z; z++) {
        int val= 0;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%d", &val);
        oField[x][y][z]= val;
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", X, Y, Z);
}


void SrtFileInput::LoadScalarFieldImagePNGFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<double>>>& oField,
    bool const iVerbose) {
  // Open the file
  std::ifstream inputFile(iFullpath.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (!inputFile.is_open()) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the file size
  std::streamsize size= 0;
  if (inputFile.seekg(0, std::ios::end).good()) size= inputFile.tellg();
  if (inputFile.seekg(0, std::ios::beg).good()) size-= inputFile.tellg();
  if (size <= 0) {
    printf("[ERROR] empty file\n\n");
    throw 0;
  }

  if (iVerbose)
    printf("Loading PNG scalar field file [%s]\n", iFullpath.c_str());

  // Load the contents of the PNG file into a buffer
  std::vector<unsigned char> buffer((size_t)size);
  inputFile.read((char*)(&buffer[0]), size);

  // Decode the buffer and store into an image vector (row/width major, col/height minor, with the 4 RGBA components for each pixel)
  unsigned long w, h;
  std::vector<unsigned char> imageVector;
  int error= SrtGenericPNGParser::DecodePNG(imageVector, w, h, buffer.empty() ? 0 : &buffer[0], (unsigned long)buffer.size(), true);
  if (error != 0) {
    printf("[ERROR] Error code %d when decoding the PNG file in buffer\n\n", error);
    throw 0;
  }

  // Convert the imageVector into the flat 3D field
  int nbY= w, nbZ= h;
  oField= std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(nbY, std::vector<double>(nbZ, 0.0)));
  for (int y= 0; y < nbY; y++) {
    for (int z= 0; z < nbZ; z++) {
      double r= double(imageVector[z * nbY * 4 + y * 4 + 0]) / 255.0;
      double g= double(imageVector[z * nbY * 4 + y * 4 + 1]) / 255.0;
      double b= double(imageVector[z * nbY * 4 + y * 4 + 2]) / 255.0;
      oField[0][y][nbZ - 1 - z]= (r + g + b) / 3.0;
    }
  }

  if (iVerbose)
    printf("File loaded: %d x %d x %d voxels\n", (int)oField.size(), (int)oField[0].size(), (int)oField[0][0].size());
}


void SrtFileInput::LoadRGBAFieldImagePNGFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<double>>>& oRField,
    std::vector<std::vector<std::vector<double>>>& oGField,
    std::vector<std::vector<std::vector<double>>>& oBField,
    std::vector<std::vector<std::vector<double>>>& oAField,
    bool const iVerbose) {
  // Open the file
  std::ifstream inputFile(iFullpath.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (!inputFile.is_open()) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the file size
  std::streamsize size= 0;
  if (inputFile.seekg(0, std::ios::end).good()) size= inputFile.tellg();
  if (inputFile.seekg(0, std::ios::beg).good()) size-= inputFile.tellg();
  if (size <= 0) {
    printf("[ERROR] empty file\n\n");
    throw 0;
  }

  if (iVerbose) {
    printf("Loading PNG scalar field file [%s]...", iFullpath.c_str());
    fflush(stdout);
  }

  // Load the contents of the PNG file into a buffer
  std::vector<unsigned char> buffer((size_t)size);
  inputFile.read((char*)(&buffer[0]), size);

  // Decode the buffer and store into an image vector (row/width major, col/height minor, with the 4 RGBA components for each pixel)
  unsigned long w, h;
  std::vector<unsigned char> imageVector;
  int error= SrtGenericPNGParser::DecodePNG(imageVector, w, h, buffer.empty() ? 0 : &buffer[0], (unsigned long)buffer.size(), true);
  if (error != 0) {
    printf("[ERROR] Error code %d when decoding the PNG file in buffer\n\n", error);
    throw 0;
  }

  // Convert the imageVector into the flat 3D field
  int nbY= w, nbZ= h;
  oRField= std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(nbY, std::vector<double>(nbZ, 0.0)));
  oGField= std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(nbY, std::vector<double>(nbZ, 0.0)));
  oBField= std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(nbY, std::vector<double>(nbZ, 0.0)));
  oAField= std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(nbY, std::vector<double>(nbZ, 0.0)));
  for (int y= 0; y < nbY; y++) {
    for (int z= 0; z < nbZ; z++) {
      oRField[0][y][nbZ - 1 - z]= double(imageVector[z * nbY * 4 + y * 4 + 0]) / 255.0;
      oGField[0][y][nbZ - 1 - z]= double(imageVector[z * nbY * 4 + y * 4 + 1]) / 255.0;
      oBField[0][y][nbZ - 1 - z]= double(imageVector[z * nbY * 4 + y * 4 + 2]) / 255.0;
      oAField[0][y][nbZ - 1 - z]= double(imageVector[z * nbY * 4 + y * 4 + 3]) / 255.0;
    }
  }

  if (iVerbose) {
    printf("File loaded: %d x %d x %d voxels\n", 1, nbY, nbZ);
  }
}


void SrtFileInput::LoadVectorFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<Eigen::Vector3d>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT vector field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int X= 0, Y= 0, Z= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &X, &Y, &Z);
  }
  if (X <= 0 || Y <= 0 || Z <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<Eigen::Vector3d>>>(X,
                                                                 std::vector<std::vector<Eigen::Vector3d>>(Y,
                                                                                                           std::vector<Eigen::Vector3d>(Z,
                                                                                                                                        Eigen::Vector3d(0.0, 0.0, 0.0))));

  // Load the field values
  for (int x= 0; x < X; x++) {
    for (int y= 0; y < Y; y++) {
      for (int z= 0; z < Z; z++) {
        double displacementX= 0.0;
        double displacementY= 0.0;
        double displacementZ= 0.0;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%lf %lf %lf", &displacementX, &displacementY, &displacementZ);
        oField[x][y][z]= {displacementX, displacementY, displacementZ};
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", X, Y, Z);
}


void SrtFileInput::LoadVectorFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<std::array<bool, 3>>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT vector field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int X= 0, Y= 0, Z= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &X, &Y, &Z);
  }
  if (X <= 0 || Y <= 0 || Z <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<std::array<bool, 3>>>>(X,
                                                                     std::vector<std::vector<std::array<bool, 3>>>(Y,
                                                                                                                   std::vector<std::array<bool, 3>>(Z,
                                                                                                                                                    {false, false, false})));

  // Load the field values
  for (int x= 0; x < X; x++) {
    for (int y= 0; y < Y; y++) {
      for (int z= 0; z < Z; z++) {
        std::array<int, 3> val;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%d %d %d", &val[0], &val[1], &val[2]);
        oField[x][y][z]= {val[0] != 0, val[1] != 0, val[2] != 0};
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", X, Y, Z);
}


void SrtFileInput::LoadVectorFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<std::array<double, 3>>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT vector field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int X= 0, Y= 0, Z= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &X, &Y, &Z);
  }
  if (X <= 0 || Y <= 0 || Z <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<std::array<double, 3>>>>(X,
                                                                       std::vector<std::vector<std::array<double, 3>>>(Y,
                                                                                                                       std::vector<std::array<double, 3>>(Z,
                                                                                                                                                          {0.0, 0.0, 0.0})));

  // Load the field values
  for (int x= 0; x < X; x++) {
    for (int y= 0; y < Y; y++) {
      for (int z= 0; z < Z; z++) {
        std::array<double, 3> val;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%lf %lf %lf", &val[0], &val[1], &val[2]);
        oField[x][y][z]= val;
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", X, Y, Z);
}


void SrtFileInput::LoadTensorFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<std::array<double, 9>>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT tensor field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int nbX= 0, nbY= 0, nbZ= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &nbX, &nbY, &nbZ);
  }
  if (nbX <= 0 || nbY <= 0 || nbZ <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<std::array<double, 9>>>>(nbX,
                                                                       std::vector<std::vector<std::array<double, 9>>>(nbY,
                                                                                                                       std::vector<std::array<double, 9>>(nbZ,
                                                                                                                                                          {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0})));

  // Load the field values
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        std::array<double, 9> val;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &val[0], &val[1], &val[2], &val[3], &val[4], &val[5], &val[6], &val[7], &val[8]);
        oField[x][y][z]= val;
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", nbX, nbY, nbZ);
}


void SrtFileInput::LoadTensorFieldTXTFile(
    std::string const iFullpath,
    std::vector<std::vector<std::vector<Eigen::Matrix3d>>>& oField,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading TXT tensor field file [%s]\n", iFullpath.c_str());

  // Open the file
  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  // Get the field dimensions
  char buffer[1000];
  int nbX= 0, nbY= 0, nbZ= 0;
  if (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    sscanf(buffer, "%d %d %d", &nbX, &nbY, &nbZ);
  }
  if (nbX <= 0 || nbY <= 0 || nbZ <= 0) {
    printf("[ERROR] Unable to read the field dimensions\n");
    throw 0;
  }

  // Allocate the field
  oField= std::vector<std::vector<std::vector<Eigen::Matrix3d>>>(nbX,
                                                                 std::vector<std::vector<Eigen::Matrix3d>>(nbY,
                                                                                                           std::vector<Eigen::Matrix3d>(nbZ)));

  // Load the field values
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        Eigen::Matrix3d val;
        if (fgets(buffer, sizeof buffer, inputFile) != NULL)
          sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &val(0), &val(1), &val(2), &val(3), &val(4), &val(5), &val(6), &val(7), &val(8));
        oField[x][y][z]= val;
      }
    }
  }

  // Close the file
  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d x %d x %d\n", nbX, nbY, nbZ);
}


void SrtFileInput::LoadMeshOBJFile(
    std::string const iFullpath,
    std::vector<Eigen::Vector3d>& oPoints,
    std::vector<Eigen::Vector3d>& oColors,
    std::vector<Eigen::Vector3i>& oTriangles,
    bool const iVerbose) {
  if (iVerbose)
    printf("Loading OBJ mesh file [%s]\n", iFullpath.c_str());

  FILE* inputFile= nullptr;
  inputFile= fopen(iFullpath.c_str(), "r");
  if (inputFile == nullptr) {
    printf("[ERROR] Unable to open the file\n\n");
    throw 0;
  }

  oPoints.clear();
  oColors.clear();
  oTriangles.clear();

  char buffer[1000];
  while (fgets(buffer, sizeof buffer, inputFile) != NULL) {
    double x, y, z;
    double r, g, b;
    int v0, v1, v2, v3;
    int n0, n1, n2, n3;
    if (sscanf(buffer, "v %lf %lf %lf %lf %lf %lf", &x, &y, &z, &r, &g, &b) == 6) {
      oPoints.push_back(Eigen::Vector3d(x, y, z));
      oColors.push_back(Eigen::Vector3d(r, g, b));
    }
    else if (sscanf(buffer, "v %lf %lf %lf", &x, &y, &z) == 3) {
      oPoints.push_back(Eigen::Vector3d(x, y, z));
      oColors.push_back(Eigen::Vector3d(0.5, 0.5, 0.5));
    }
    else if (sscanf(buffer, "f %d//%d %d//%d %d//%d %d//%d", &v0, &n0, &v1, &n1, &v2, &n2, &v3, &n3) == 8) {
      oTriangles.push_back(Eigen::Vector3i(v0 - 1, v1 - 1, v2 - 1));
      oTriangles.push_back(Eigen::Vector3i(v0 - 1, v2 - 1, v3 - 1));
    }
    else if (sscanf(buffer, "f %d//%d %d//%d %d//%d", &v0, &n0, &v1, &n1, &v2, &n2) == 6) {
      oTriangles.push_back(Eigen::Vector3i(v0 - 1, v1 - 1, v2 - 1));
    }
    else if (sscanf(buffer, "f %d %d %d %d", &v0, &v1, &v2, &v3) == 4) {
      oTriangles.push_back(Eigen::Vector3i(v0 - 1, v1 - 1, v2 - 1));
      oTriangles.push_back(Eigen::Vector3i(v0 - 1, v2 - 1, v3 - 1));
    }
    else if (sscanf(buffer, "f %d %d %d", &v0, &v1, &v2) == 3) {
      oTriangles.push_back(Eigen::Vector3i(v0 - 1, v1 - 1, v2 - 1));
    }
  }

  fclose(inputFile);

  if (iVerbose)
    printf("File loaded: %d points, %d triangles\n", int(oPoints.size()), int(oTriangles.size()));
}
