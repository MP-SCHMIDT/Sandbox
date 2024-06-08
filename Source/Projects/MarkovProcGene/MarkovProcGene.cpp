#include "MarkovProcGene.hpp"


// Standard lib
#include <array>
#include <cmath>
#include <vector>

// GLUT lib
#include "GL/freeglut.h"

// Algo headers
#include "Draw/DrawShape.hpp"
#include "FileIO/FileInput.hpp"
#include "Math/Field.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
MarkovProcGene::MarkovProcGene() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void MarkovProcGene::SetActiveProject() {
  if (!isActivProj || D.UI.empty()) {
    D.UI.clear();
    D.UI.push_back(ParamUI("Scenario________", 0));
    D.UI.push_back(ParamUI("ResolutionX_____", 1));
    D.UI.push_back(ParamUI("ResolutionY_____", 21));
    D.UI.push_back(ParamUI("ResolutionZ_____", 21));
    D.UI.push_back(ParamUI("RuleSizeX_______", 4));
    D.UI.push_back(ParamUI("RuleSizeY_______", 4));
    D.UI.push_back(ParamUI("RuleSizeZ_______", 4));
    D.UI.push_back(ParamUI("NbSubsti________", 1));
    D.UI.push_back(ParamUI("ShadeCoeff______", 1));
    D.UI.push_back(ParamUI("VerboseLevel____", 0));
  }

  if (D.UI.size() != VerboseLevel____ + 1) {
    printf("[ERROR] Invalid parameter count in UI\n");
  }

  D.boxMin= {0.0, 0.0, 0.0};
  D.boxMax= {1.0, 1.0, 1.0};

  isActivProj= true;
  isAllocated= false;
  isRefreshed= false;
}


// Check if parameter changes should trigger an allocation
bool MarkovProcGene::CheckAlloc() {
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool MarkovProcGene::CheckRefresh() {
  if (D.UI[Scenario________].hasChanged()) isRefreshed= false;
  if (D.UI[ResolutionX_____].hasChanged()) isRefreshed= false;
  if (D.UI[ResolutionY_____].hasChanged()) isRefreshed= false;
  if (D.UI[ResolutionZ_____].hasChanged()) isRefreshed= false;
  if (D.UI[RuleSizeX_______].hasChanged()) isRefreshed= false;
  if (D.UI[RuleSizeY_______].hasChanged()) isRefreshed= false;
  if (D.UI[RuleSizeZ_______].hasChanged()) isRefreshed= false;
  return isRefreshed;
}


// Allocate the project data
void MarkovProcGene::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;
}


// Refresh the project
void MarkovProcGene::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Get UI parameters
  nbX= std::max(D.UI[ResolutionX_____].I(), 1);
  nbY= std::max(D.UI[ResolutionY_____].I(), 1);
  nbZ= std::max(D.UI[ResolutionZ_____].I(), 1);

  // Reset progress counters
  activeSet= -1;
  activeRul= -1;
  D.Plot.clear();

  // Initialize dictionnary and field values
  Dict.clear();
  Field= Field::AllocField3D(nbX, nbY, nbZ, 0);
  std::array<std::vector<std::vector<std::vector<int>>>, 2> tmpRule;
  int scenario= 0;

  // Wave Function collapse from BMP image
  if (scenario++ == D.UI[Scenario________].I()) {
    nbX= 2;
    Field= Field::AllocField3D(nbX, nbY, nbZ, 0);

    std::vector<std::vector<std::array<float, 4>>> imageRGBA;
    FileInput::LoadImageBMPFile("./FileInput/Images/WFC_Example.bmp", imageRGBA, false);
    const int nbWImag= (int)imageRGBA.size();
    const int nbHImag= (int)imageRGBA[0].size();
    const int nbWCell= std::min(std::max(D.UI[RuleSizeY_______].I(), 3), nbWImag);
    const int nbHCell= std::min(std::max(D.UI[RuleSizeZ_______].I(), 3), nbHImag);

    std::vector<std::vector<int>> Imag= Field::AllocField2D(nbWImag, nbHImag, 0);
    for (int wImag= 0; wImag < nbWImag; wImag++) {
      for (int hImag= 0; hImag < nbHImag; hImag++) {
        if (imageRGBA[wImag][hImag][3] < 0.5f) Imag[wImag][hImag]= 0;
        else if (imageRGBA[wImag][hImag][0] > 0.5f && imageRGBA[wImag][hImag][1] < 0.5f && imageRGBA[wImag][hImag][2] < 0.5f) Imag[wImag][hImag]= 1;
        else if (imageRGBA[wImag][hImag][0] < 0.5f && imageRGBA[wImag][hImag][1] > 0.5f && imageRGBA[wImag][hImag][2] < 0.5f) Imag[wImag][hImag]= 2;
        else if (imageRGBA[wImag][hImag][0] < 0.5f && imageRGBA[wImag][hImag][1] < 0.5f && imageRGBA[wImag][hImag][2] > 0.5f) Imag[wImag][hImag]= 3;
        else if (imageRGBA[wImag][hImag][0] < 0.5f && imageRGBA[wImag][hImag][1] > 0.5f && imageRGBA[wImag][hImag][2] > 0.5f) Imag[wImag][hImag]= 4;
        else if (imageRGBA[wImag][hImag][0] > 0.5f && imageRGBA[wImag][hImag][1] < 0.5f && imageRGBA[wImag][hImag][2] > 0.5f) Imag[wImag][hImag]= 5;
        else if (imageRGBA[wImag][hImag][0] > 0.5f && imageRGBA[wImag][hImag][1] > 0.5f && imageRGBA[wImag][hImag][2] < 0.5f) Imag[wImag][hImag]= 6;
        else if (imageRGBA[wImag][hImag][0] < 0.5f && imageRGBA[wImag][hImag][1] < 0.5f && imageRGBA[wImag][hImag][2] < 0.5f) Imag[wImag][hImag]= 7;
        else if (imageRGBA[wImag][hImag][0] > 0.5f && imageRGBA[wImag][hImag][1] > 0.5f && imageRGBA[wImag][hImag][2] > 0.5f) Imag[wImag][hImag]= 9;
        else Imag[wImag][hImag]= 8;
      }
    }

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    for (int wImag= 0; wImag < nbWImag - (nbWCell - 1); wImag++) {
      for (int hImag= 0; hImag < nbHImag - (nbHCell - 1); hImag++) {
        // for (int wImag= 0; wImag < nbWImag - (nbWCell - 1); wImag+= nbWCell - 1) {
        //   for (int hImag= 0; hImag < nbHImag - (nbHCell - 1); hImag+= nbHCell - 1) {
        for (int useT= 0; useT < 2; useT++) {
          for (int useB= 0; useB < 2; useB++) {
            for (int useR= 0; useR < 2; useR++) {
              for (int useL= 0; useL < 2; useL++) {
                if (useT == 0 && useB == 0 && useR == 0 && useL == 0) continue;
                tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(2, nbWCell, nbHCell, 0), Field::AllocField3D(2, nbWCell, nbHCell, 0)});
                if (useB > 0 || useL > 0) tmpRule[0][0][0][0]= 1;
                if (useT > 0 || useL > 0) tmpRule[0][0][0][nbHCell - 1]= 1;
                if (useB > 0 || useR > 0) tmpRule[0][0][nbWCell - 1][0]= 1;
                if (useT > 0 || useR > 0) tmpRule[0][0][nbWCell - 1][nbHCell - 1]= 1;
                tmpRule[1][0][0][0]= 1;
                tmpRule[1][0][0][nbHCell - 1]= 1;
                tmpRule[1][0][nbWCell - 1][0]= 1;
                tmpRule[1][0][nbWCell - 1][nbHCell - 1]= 1;
                for (int idxW= 0; idxW < nbWCell; idxW++) {
                  for (int idxH= 0; idxH < nbHCell; idxH++) {
                    if (useT > 0 && idxH == nbHCell - 1) tmpRule[0][1][idxW][idxH]= Imag[wImag + idxW][hImag + idxH];
                    if (useB > 0 && idxH == 0) tmpRule[0][1][idxW][idxH]= Imag[wImag + idxW][hImag + idxH];
                    if (useR > 0 && idxW == nbWCell - 1) tmpRule[0][1][idxW][idxH]= Imag[wImag + idxW][hImag + idxH];
                    if (useL > 0 && idxW == 0) tmpRule[0][1][idxW][idxH]= Imag[wImag + idxW][hImag + idxH];
                    tmpRule[1][1][idxW][idxH]= Imag[wImag + idxW][hImag + idxH];
                  }
                }

                Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
                // Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, -3, tmpRule));
                // Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, -2, +3, tmpRule));
                // Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, -2, -3, tmpRule));
              }
            }
          }
        }
      }
    }

    if (nbWCell <= nbY && nbHCell <= nbZ) {
      int idxSeedCellW= Random::Val(0, (nbY - 1) - (nbWCell - 1));
      int idxSeedCellH= Random::Val(0, (nbZ - 1) - (nbHCell - 1));
      int idxSeedRule= Random::Val(0, (int)Dict[0].size() - 1);
      for (int x= 0; x < 2; x++)
        for (int idxW= 0; idxW < (int)Dict[0][idxSeedRule][0][0].size(); idxW++)
          for (int idxH= 0; idxH < (int)Dict[0][idxSeedRule][0][0][0].size(); idxH++)
            Field[x][idxSeedCellW + idxW][idxSeedCellH + idxH]= Dict[0][idxSeedRule][1][x][idxW][idxH];
    }
  }

  // Random noise
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 1, 0), Field::AllocField3D(1, 1, 1, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 0;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 1;
  }

  // Random noise with multiple sets
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 1, 0), Field::AllocField3D(1, 1, 1, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 0;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 1;

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 1, 0), Field::AllocField3D(1, 1, 1, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 1;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 2;
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 1, 0), Field::AllocField3D(1, 1, 1, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 1;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 3;

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 1, 0), Field::AllocField3D(1, 1, 1, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 2;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 1;
  }

  // Spread
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 2, 0), Field::AllocField3D(1, 1, 2, 0)});
    tmpRule[0][0][0][0]= 0;
    tmpRule[0][0][0][1]= 1;
    tmpRule[1][0][0][0]= 1;
    tmpRule[1][0][0][1]= 1;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, -3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, -3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+3, +1, +2, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-3, +1, +2, tmpRule));

    Field[nbX / 2][nbY / 2][nbZ / 2]= 1;
  }

  // Tron infinite wall
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 3, 0), Field::AllocField3D(1, 1, 3, 0)});
    tmpRule[0][0][0][0]= 1;
    tmpRule[1][0][0][0]= 2;
    tmpRule[1][0][0][1]= 2;
    tmpRule[1][0][0][2]= 1;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, -3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, -3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+3, +1, +2, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-3, +1, +2, tmpRule));

    Field[nbX / 2][nbY / 2][nbZ / 2]= 1;
  }

  // Snake
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 5, 0), Field::AllocField3D(1, 1, 5, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][2]= 3;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 4;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][1]= 2;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][2]= 2;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][3]= 2;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][4]= 1;

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 3, 0), Field::AllocField3D(1, 1, 3, 0)});
    tmpRule[0][0][0][0]= 1;
    tmpRule[1][0][0][0]= 2;
    tmpRule[1][0][0][1]= 2;
    tmpRule[1][0][0][2]= 1;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, -3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, -3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+3, +1, +2, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-3, +1, +2, tmpRule));

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 3, 0), Field::AllocField3D(1, 1, 3, 0)});
    tmpRule[0][0][0][0]= 4;
    tmpRule[0][0][0][1]= 2;
    tmpRule[0][0][0][2]= 2;
    tmpRule[1][0][0][2]= 4;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, -3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, -3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+3, +1, +2, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-3, +1, +2, tmpRule));

    Field[nbX / 2][nbY / 2][nbZ / 2]= 3;
  }

  // Spanning tree
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 3, 0), Field::AllocField3D(1, 1, 3, 0)});
    tmpRule[0][0][0][0]= 3;
    tmpRule[1][0][0][0]= 3;
    tmpRule[1][0][0][1]= 8;
    tmpRule[1][0][0][2]= 3;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, -3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, -3, +1, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+3, +1, +2, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-3, +1, +2, tmpRule));

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 1, 0), Field::AllocField3D(1, 1, 1, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 3;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 8;

    Field[nbX / 2][nbY / 2][nbZ / 2]= 3;
  }

  // Flower garden
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(5, 5, 1, 0), Field::AllocField3D(5, 5, 1, 0)}));
    FillRuleBox(Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1], 0, 0, 0, 4, 4, 0, 2, true, true);
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][2][2][0]= 2;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][2][2][0]= 1;

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 2, 0), Field::AllocField3D(1, 1, 2, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][0]= 1;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 1;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][1]= 6;

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(5, 5, 2, 0), Field::AllocField3D(5, 5, 2, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][2][2][0]= 6;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][2][2][0]= 2;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][2][2][1]= 6;

    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(5, 5, 2, 0), Field::AllocField3D(5, 5, 2, 0)});
    tmpRule[0][2][2][0]= 6;
    tmpRule[1][2][2][0]= 2;
    tmpRule[1][2][1][1]= 6;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, -2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +1, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-2, +1, +3, tmpRule));

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(3, 3, 1, 0), Field::AllocField3D(3, 3, 1, 0)});
    tmpRule[0][1][1][0]= 6;
    tmpRule[1][1][1][0]= 6;
    tmpRule[1][1][0][0]= 7;
    tmpRule[1][0][1][0]= 7;
    tmpRule[1][1][2][0]= 7;
    tmpRule[1][2][1][0]= 7;
    Dict[(int)Dict.size() - 1].push_back(BuildColorSwap(7, 3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildColorSwap(7, 4, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildColorSwap(7, 5, tmpRule));

    for (int x= 0; x < nbX; x++)
      for (int y= 0; y < nbY; y++)
        Field[x][y][0]= 2;
  }

  // Galton board
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 2, 0), Field::AllocField3D(1, 1, 2, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][1]= 3;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][1]= 3;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 6;
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 2, 2, 0), Field::AllocField3D(1, 2, 2, 0)});
    tmpRule[0][0][0][1]= 4;
    tmpRule[0][0][1][0]= 4;
    tmpRule[0][0][1][1]= 6;
    tmpRule[1][0][0][0]= 6;
    tmpRule[1][0][0][1]= 4;
    tmpRule[1][0][1][0]= 4;
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, -2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +1, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-2, +1, +3, tmpRule));

    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 2, 0), Field::AllocField3D(1, 1, 2, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][1]= 6;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 6;
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 3, 0), Field::AllocField3D(1, 1, 3, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][2]= 6;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 6;
    Dict[(int)Dict.size() - 1].push_back(std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 1, 4, 0), Field::AllocField3D(1, 1, 4, 0)}));
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][0][0][0][3]= 6;
    Dict[(int)Dict.size() - 1][(int)Dict[Dict.size() - 1].size() - 1][1][0][0][0]= 6;

    for (int x= 0; x < nbX; x++) {
      for (int y= 0; y < nbY; y++) {
        for (int z= 0; z < nbZ; z++) {
          if ((x + y + z) % 2 == 0 && z > 3 * nbZ / 5 && z < nbZ - 2)
            Field[x][y][z]= 4;
          if (std::abs(x - (nbX - 1) / 2) <= 1 && std::abs(y - (nbY - 1) / 2) <= 1 && z == nbZ - 1)
            Field[x][y][z]= 3;
        }
      }
    }
  }

  // Skyscraper
  if (scenario++ == D.UI[Scenario________].I()) {
    int const blocX= std::max(D.UI[RuleSizeX_______].I(), 3);
    int const blocY= blocX;
    int const blocZ= std::max(D.UI[RuleSizeZ_______].I(), 3);

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(blocX, blocY, blocZ, 0), Field::AllocField3D(blocX, blocY, blocZ, 0)});
    FillRuleBox(tmpRule, 0, 0, 0, blocX - 1, blocY - 1, 0, 8, true, true);
    FillRuleBox(tmpRule, 0, 0, blocZ - 1, blocX - 1, blocY - 1, blocZ - 1, 8, false, true);

    FillRuleBox(tmpRule, 0, 0, blocZ - 1, 0, blocY - 1, blocZ - 1, 8, true, true);

    FillRuleBox(tmpRule, 0, 0, 1, 0, 0, blocZ - 2, 8, true, true);
    FillRuleBox(tmpRule, 0, blocY - 1, 1, 0, blocY - 1, blocZ - 2, 8, true, true);
    FillRuleBox(tmpRule, blocX - 1, 0, 1, blocX - 1, 0, blocZ - 2, 8, false, true);
    FillRuleBox(tmpRule, blocX - 1, blocY - 1, 1, blocX - 1, blocY - 1, blocZ - 2, 8, false, true);
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(-1, +2, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, +1, +3, tmpRule));
    Dict[(int)Dict.size() - 1].push_back(BuildSymmetric(+2, -1, +3, tmpRule));

    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(blocX, blocY, blocZ, 0), Field::AllocField3D(blocX, blocY, blocZ, 0)});
    FillRuleBox(tmpRule, 0, 0, 0, blocX - 1, blocY - 1, 0, 8, true, true);
    FillRuleBox(tmpRule, 0, 0, blocZ - 1, blocX - 1, blocY - 1, blocZ - 1, 8, false, true);

    FillRuleBox(tmpRule, 0, 0, 1, 0, 0, blocZ - 2, 8, false, true);
    FillRuleBox(tmpRule, 0, blocY - 1, 1, 0, blocY - 1, blocZ - 2, 8, false, true);
    FillRuleBox(tmpRule, blocX - 1, 0, 1, blocX - 1, 0, blocZ - 2, 8, false, true);
    FillRuleBox(tmpRule, blocX - 1, blocY - 1, 1, blocX - 1, blocY - 1, blocZ - 2, 8, false, true);

    Dict[(int)Dict.size() - 1].push_back(tmpRule);
    FillRuleBox(tmpRule, 1, 1, blocZ - 1, blocX - 2, blocY - 2, blocZ - 1, 0, false, true);
    Dict[(int)Dict.size() - 1].push_back(tmpRule);

    for (int x= 0; x < nbX; x++)
      for (int y= 0; y < nbY; y++)
        Field[x][y][0]= 8;
  }


  // Conway game of life
  if (scenario++ == D.UI[Scenario________].I()) {
    Dict.push_back(std::vector<std::array<std::vector<std::vector<std::vector<int>>>, 2>>());
    for (int useCC= 0; useCC < 2; useCC++) {
      for (int useNN= 0; useNN < 2; useNN++) {
        for (int useNE= 0; useNE < 2; useNE++) {
          for (int useEE= 0; useEE < 2; useEE++) {
            for (int useSE= 0; useSE < 2; useSE++) {
              for (int useSS= 0; useSS < 2; useSS++) {
                for (int useSO= 0; useSO < 2; useSO++) {
                  for (int useOO= 0; useOO < 2; useOO++) {
                    for (int useNO= 0; useNO < 2; useNO++) {
                      tmpRule= std::array<std::vector<std::vector<std::vector<int>>>, 2>({Field::AllocField3D(1, 3, 3, 0), Field::AllocField3D(1, 3, 3, 0)});
                      if (useCC == 1) tmpRule[0][0][1][1]= 8;
                      if (useNN == 1) tmpRule[0][0][1][2]= 8;
                      if (useNE == 1) tmpRule[0][0][2][2]= 8;
                      if (useEE == 1) tmpRule[0][0][2][1]= 8;
                      if (useSE == 1) tmpRule[0][0][2][0]= 8;
                      if (useSS == 1) tmpRule[0][0][1][0]= 8;
                      if (useSO == 1) tmpRule[0][0][0][0]= 8;
                      if (useOO == 1) tmpRule[0][0][0][1]= 8;
                      if (useNO == 1) tmpRule[0][0][0][2]= 8;
                      tmpRule[1]= tmpRule[0];
                      if (useCC == 1 && useNN + useNE + useEE + useSE + useSS + useSO + useOO + useNO < 2) {
                        tmpRule[1][0][1][1]= 0;
                        Dict[(int)Dict.size() - 1].push_back(tmpRule);
                      }
                      else if (useCC == 1 && useNN + useNE + useEE + useSE + useSS + useSO + useOO + useNO > 3) {
                        tmpRule[1][0][1][1]= 0;
                        Dict[(int)Dict.size() - 1].push_back(tmpRule);
                      }
                      else if (useCC == 0 && useNN + useNE + useEE + useSE + useSS + useSO + useOO + useNO == 3) {
                        tmpRule[1][0][1][1]= 8;
                        Dict[(int)Dict.size() - 1].push_back(tmpRule);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    for (int x= 0; x < nbX; x++)
      for (int y= 0; y < nbY; y++)
        for (int z= 0; z < nbZ; z++)
          if (Random::Val(0, 3) == 0)
            Field[x][y][z]= 8;
  }


  // Check and fix for negative signs in field or rules
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        if (Field[x][y][z] < 0) {
          printf("fixed negative fuckup in field\n");
          Field[x][y][z]= std::abs(Field[x][y][z]);
        }
      }
    }
  }
  for (int idxSet= 0; idxSet < (int)Dict.size(); idxSet++) {
    for (int idxRule= 0; idxRule < (int)Dict[idxSet].size(); idxRule++) {
      for (int k= 0; k < 2; k++) {
        for (int xR= 0; xR < (int)Dict[idxSet][idxRule][0].size(); xR++) {
          for (int yR= 0; yR < (int)Dict[idxSet][idxRule][0][0].size(); yR++) {
            for (int zR= 0; zR < (int)Dict[idxSet][idxRule][0][0][0].size(); zR++) {
              if (Dict[idxSet][idxRule][k][xR][yR][zR] < 0) {
                printf("fixed negative fuckup in rules\n");
                Dict[idxSet][idxRule][k][xR][yR][zR]= std::abs(Dict[idxSet][idxRule][k][xR][yR][zR]);
              }
            }
          }
        }
      }
    }
  }
}


// Handle keypress
void MarkovProcGene::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Handle mouse action
void MarkovProcGene::MousePress(const unsigned char mouse) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)mouse;  // Disable warning unused variable
}


// Animate the project
void MarkovProcGene::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Iterate over the desired number of substitutions
  if (Dict.empty()) return;
  for (int idxIter= 0; idxIter < D.UI[NbSubsti________].I(); idxIter++) {
    activeRul= -1;
    activeSet= -1;
    int activeMatchCount= 0;

    for (int idxSet= 0; idxSet < (int)Dict.size(); idxSet++) {
      int matchCount= 0;
      for (int idxRule= 0; idxRule < (int)Dict[idxSet].size(); idxRule++) {
        int nbXRule= (int)Dict[idxSet][idxRule][0].size();
        int nbYRule= (int)Dict[idxSet][idxRule][0][0].size();
        int nbZRule= (int)Dict[idxSet][idxRule][0][0][0].size();
        for (int xF= 0; xF <= nbX - nbXRule; xF++) {
          for (int yF= 0; yF <= nbY - nbYRule; yF++) {
            for (int zF= 0; zF <= nbZ - nbZRule; zF++) {
              bool isMatch= true;
              for (int xR= 0; xR < nbXRule && isMatch; xR++)
                for (int yR= 0; yR < nbYRule && isMatch; yR++)
                  for (int zR= 0; zR < nbZRule && isMatch; zR++)
                    // if (Field[xF + xR][yF + yR][zF + zR] != 0)
                    if (Field[xF + xR][yF + yR][zF + zR] != Dict[idxSet][idxRule][0][xR][yR][zR])
                      isMatch= false;
              if (isMatch)
                matchCount++;
            }
          }
        }
      }
      if (activeSet < 0 && matchCount > 0) {
        activeSet= idxSet;
        activeMatchCount= matchCount;
      }
      D.Plot.resize(Dict.size());
      D.Plot[idxSet].name= "Set" + std::to_string(idxSet);
      D.Plot[idxSet].val.push_back((double)matchCount);
    }

    if (activeSet < 0) continue;

    // Choose a random match to substitute according to its rule
    int matchChosen= Random::Val(0, activeMatchCount);
    bool substitutionDone= false;
    for (int idxRule= 0; idxRule < (int)Dict[activeSet].size() && !substitutionDone; idxRule++) {
      int nbXRule= (int)Dict[activeSet][idxRule][0].size();
      int nbYRule= (int)Dict[activeSet][idxRule][0][0].size();
      int nbZRule= (int)Dict[activeSet][idxRule][0][0][0].size();
      for (int xF= 0; xF <= nbX - nbXRule && !substitutionDone; xF++) {
        for (int yF= 0; yF <= nbY - nbYRule && !substitutionDone; yF++) {
          for (int zF= 0; zF <= nbZ - nbZRule && !substitutionDone; zF++) {
            bool isMatch= true;
            for (int xR= 0; xR < nbXRule && isMatch; xR++)
              for (int yR= 0; yR < nbYRule && isMatch; yR++)
                for (int zR= 0; zR < nbZRule && isMatch; zR++)
                  // if (Field[xF + xR][yF + yR][zF + zR] != 0)
                  if (Field[xF + xR][yF + yR][zF + zR] != Dict[activeSet][idxRule][0][xR][yR][zR])
                    isMatch= false;
            if (isMatch)
              activeMatchCount--;
            if (activeMatchCount == matchChosen) {
              for (int xR= 0; xR < nbXRule && isMatch; xR++)
                for (int yR= 0; yR < nbYRule && isMatch; yR++)
                  for (int zR= 0; zR < nbZRule && isMatch; zR++)
                    Field[xF + xR][yF + yR][zF + zR]= Dict[activeSet][idxRule][1][xR][yR][zR];
              substitutionDone= true;
              activeRul= idxRule;
            }
          }
        }
      }
    }
  }
}


void util_SetColorVoxel(const int iVal, const float iShading) {
  if (iVal == 0) glColor3f(0.0f, 0.0f, 0.0f);
  if (iVal == 1) glColor3f(iShading * 0.8f, iShading * 0.5f, iShading * 0.5f);
  if (iVal == 2) glColor3f(iShading * 0.5f, iShading * 0.8f, iShading * 0.5f);
  if (iVal == 3) glColor3f(iShading * 0.5f, iShading * 0.5f, iShading * 0.8f);
  if (iVal == 4) glColor3f(iShading * 0.3f, iShading * 0.7f, iShading * 0.7f);
  if (iVal == 5) glColor3f(iShading * 0.7f, iShading * 0.3f, iShading * 0.7f);
  if (iVal == 6) glColor3f(iShading * 0.7f, iShading * 0.7f, iShading * 0.3f);
  if (iVal == 7) glColor3f(iShading * 0.3f, iShading * 0.3f, iShading * 0.3f);
  if (iVal == 8) glColor3f(iShading * 0.7f, iShading * 0.7f, iShading * 0.7f);
  if (iVal == 9) glColor3f(iShading * 1.0f, iShading * 1.0f, iShading * 1.0f);
}


// Draw the project
void MarkovProcGene::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Get sizes
  int maxDim= std::max(std::max(nbX, nbY), nbZ);
  float voxSize= 1.0 / (float)maxDim;

  // Compute the shading directions
  std::vector<std::array<int, 3>> ShaDir;
  for (int x= -1; x <= 1; x++)
    for (int y= -1; y <= 1; y++)
      for (int z= -1; z <= 1; z++)
        if (x != 0 || y != 0 || z != 0)
          ShaDir.push_back(std::array<int, 3>({x, y, z}));

  // Compute the voxel shading map
  std::vector<std::vector<std::vector<float>>> FieldVisi= Field::AllocField3D(nbX, nbY, nbZ, 0.0f);
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        int countVisi= 0;
        for (int idxDir= 0; idxDir < (int)ShaDir.size(); idxDir++) {
          int xOff= x;
          int yOff= y;
          int zOff= z;
          bool isOccluded= false;
          while (!isOccluded) {
            xOff+= ShaDir[idxDir][0];
            yOff+= ShaDir[idxDir][1];
            zOff+= ShaDir[idxDir][2];
            if (xOff < 0 || xOff >= nbX) break;
            if (yOff < 0 || yOff >= nbY) break;
            if (zOff < 0 || zOff >= nbZ) break;
            if (Field[xOff][yOff][zOff] > 0) isOccluded= true;
          }
          if (!isOccluded)
            countVisi++;
        }
        FieldVisi[x][y][z]= 0.5f + 0.5f * (float)countVisi / (float)ShaDir.size();
      }
    }
  }

  // Smooth the voxel shading map
  std::vector<std::vector<std::vector<float>>> FieldVisiOld= FieldVisi;
  for (int x= 0; x < nbX; x++) {
    for (int y= 0; y < nbY; y++) {
      for (int z= 0; z < nbZ; z++) {
        FieldVisi[x][y][z]= 0.0f;
        int count= 0;
        for (int xOff= std::max(x - 1, 0); xOff <= std::min(x + 1, nbX - 1); xOff++) {
          for (int yOff= std::max(y - 1, 0); yOff <= std::min(y + 1, nbY - 1); yOff++) {
            for (int zOff= std::max(z - 1, 0); zOff <= std::min(z + 1, nbZ - 1); zOff++) {
              if (Field[xOff][yOff][zOff] > 0) {
                FieldVisi[x][y][z]+= FieldVisiOld[xOff][yOff][zOff];
                count++;
              }
            }
          }
        }
        FieldVisi[x][y][z]/= (float)count;
      }
    }
  }

  // Draw the voxels
  if (D.displayMode1) {
    glEnable(GL_LIGHTING);
    for (int x= 0; x < nbX; x++) {
      for (int y= 0; y < nbY; y++) {
        for (int z= 0; z < nbZ; z++) {
          if (Field[x][y][z] > 0) {
            util_SetColorVoxel(Field[x][y][z], 1.0f - D.UI[ShadeCoeff______].F() * (1.0f - FieldVisi[x][y][z]));
            DrawShape::DrawBoxPosSiz(0.5f - 0.5f * (float)nbX / (float)maxDim + (float)x * voxSize,
                                     0.5f - 0.5f * (float)nbY / (float)maxDim + (float)y * voxSize,
                                     0.5f - 0.5f * (float)nbZ / (float)maxDim + (float)z * voxSize,
                                     voxSize, voxSize, voxSize, true);
          }
        }
      }
    }
    glDisable(GL_LIGHTING);
  }

  // Draw the dictionnary
  if (D.displayMode2) {
    glLineWidth(3.0);
    int curOffsetY= 1;
    int currentRul= 0;
    for (int idxSet= 0; idxSet < (int)Dict.size(); idxSet++) {
      int curOffsetZ= 0;
      int maxOffsetY= 0;
      for (int idxRule= 0; idxRule < (int)Dict[idxSet].size(); idxRule++) {
        int nbXRule= (int)Dict[idxSet][idxRule][0].size();
        int nbYRule= (int)Dict[idxSet][idxRule][0][0].size();
        int nbZRule= (int)Dict[idxSet][idxRule][0][0][0].size();
        float begXI= 0.5f;
        float begYI= 1.0f + curOffsetY * voxSize;
        float begZI= 0.0f + curOffsetZ * voxSize;
        float begXO= 0.5f;
        float begYO= 1.0f + curOffsetY * voxSize + (nbYRule + 1) * voxSize;
        float begZO= 0.0f + curOffsetZ * voxSize;
        if (idxSet == activeSet && idxRule == activeRul)
          glColor3f(0.8f, 0.8f, 0.8f);
        else
          glColor3f(0.3f, 0.3f, 0.3f);
        DrawShape::DrawBoxPosSiz(begXI, begYI, begZI, nbXRule * voxSize, nbYRule * voxSize, nbZRule * voxSize, false);
        DrawShape::DrawBoxPosSiz(begXO, begYO, begZO, nbXRule * voxSize, nbYRule * voxSize, nbZRule * voxSize, false);
        glEnable(GL_LIGHTING);
        for (int xR= 0; xR < nbXRule; xR++) {
          for (int yR= 0; yR < nbYRule; yR++) {
            for (int zR= 0; zR < nbZRule; zR++) {
              if (Dict[idxSet][idxRule][0][xR][yR][zR] != 0) {
                util_SetColorVoxel(Dict[idxSet][idxRule][0][xR][yR][zR], 0.8f);
                DrawShape::DrawBoxPosSiz(begXI + xR * voxSize, begYI + yR * voxSize, begZI + zR * voxSize, voxSize, voxSize, voxSize, true);
              }
              if (Dict[idxSet][idxRule][1][xR][yR][zR] != 0) {
                util_SetColorVoxel(Dict[idxSet][idxRule][1][xR][yR][zR], 0.8f);
                DrawShape::DrawBoxPosSiz(begXO + xR * voxSize, begYO + yR * voxSize, begZO + zR * voxSize, voxSize, voxSize, voxSize, true);
              }
            }
          }
        }
        glDisable(GL_LIGHTING);
        maxOffsetY= std::max(maxOffsetY, nbYRule);
        curOffsetZ+= nbZRule + 1;
        currentRul++;
      }
      curOffsetY+= 2 * (maxOffsetY + 1) + 2;
    }
    glLineWidth(1.0);
  }

  // Draw the box
  if (D.displayMode3) {
    glLineWidth(3.0);
    glColor3f(0.5f, 0.5f, 0.5f);
    DrawShape::DrawBoxPosSiz(0.5f - 0.5f * (float)nbX / (float)maxDim,
                             0.5f - 0.5f * (float)nbY / (float)maxDim,
                             0.5f - 0.5f * (float)nbZ / (float)maxDim,
                             voxSize * nbX, voxSize * nbY, voxSize * nbZ, false);
    glLineWidth(1.0);
  }
}


void MarkovProcGene::FillRuleBox(std::array<std::vector<std::vector<std::vector<int>>>, 2>& ioRule,
                                 const int iMinX, const int iMinY, const int iMinZ,
                                 const int iMaxX, const int iMaxY, const int iMaxZ,
                                 const int iVal, const bool iFillI, const bool iFillO) {
  int nbXRule= (int)ioRule[0].size();
  int nbYRule= (int)ioRule[0][0].size();
  int nbZRule= (int)ioRule[0][0][0].size();
  for (int x= std::max(iMinX, 0); x <= std::min(iMaxX, nbXRule - 1); x++) {
    for (int y= std::max(iMinY, 0); y <= std::min(iMaxY, nbYRule - 1); y++) {
      for (int z= std::max(iMinZ, 0); z <= std::min(iMaxZ, nbZRule - 1); z++) {
        if (iFillI) ioRule[0][x][y][z]= iVal;
        if (iFillO) ioRule[1][x][y][z]= iVal;
      }
    }
  }
}


std::array<std::vector<std::vector<std::vector<int>>>, 2> MarkovProcGene::BuildSymmetric(const int iDim1, const int iDim2, const int iDim3,
                                                                                         const std::array<std::vector<std::vector<std::vector<int>>>, 2>& iRule) {
  if (iRule[0].empty()) throw;
  if (iRule[0][0].empty()) throw;
  if (iRule[0][0][0].empty()) throw;
  int nbXS= (int)iRule[0].size();
  int nbYS= (int)iRule[0][0].size();
  int nbZS= (int)iRule[0][0][0].size();
  if (iDim1 == 0 || std::abs(iDim1) > 3) throw;
  if (iDim2 == 0 || std::abs(iDim2) > 3) throw;
  if (iDim3 == 0 || std::abs(iDim3) > 3) throw;
  if (std::abs(iDim1) == std::abs(iDim2)) throw;
  if (std::abs(iDim1) == std::abs(iDim3)) throw;
  if (std::abs(iDim2) == std::abs(iDim3)) throw;

  int nbXD= (std::abs(iDim1) == 1) ? (nbXS) : ((std::abs(iDim1) == 2) ? (nbYS) : (nbZS));
  int nbYD= (std::abs(iDim2) == 1) ? (nbXS) : ((std::abs(iDim2) == 2) ? (nbYS) : (nbZS));
  int nbZD= (std::abs(iDim3) == 1) ? (nbXS) : ((std::abs(iDim3) == 2) ? (nbYS) : (nbZS));
  std::array<std::vector<std::vector<std::vector<int>>>, 2> oRule({Field::AllocField3D(nbXD, nbYD, nbZD, 0), Field::AllocField3D(nbXD, nbYD, nbZD, 0)});
  for (int xS= 0; xS < nbXS; xS++) {
    for (int yS= 0; yS < nbYS; yS++) {
      for (int zS= 0; zS < nbZS; zS++) {
        int xD= 0, yD= 0, zD= 0;
        if (iDim1 < 0) xD= (std::abs(iDim1) == 1) ? (nbXS - 1 - xS) : ((std::abs(iDim1) == 2) ? (nbYS - 1 - yS) : (nbZS - 1 - zS));
        if (iDim1 > 0) xD= (std::abs(iDim1) == 1) ? (xS) : ((std::abs(iDim1) == 2) ? (yS) : (zS));
        if (iDim2 < 0) yD= (std::abs(iDim2) == 1) ? (nbXS - 1 - xS) : ((std::abs(iDim2) == 2) ? (nbYS - 1 - yS) : (nbZS - 1 - zS));
        if (iDim2 > 0) yD= (std::abs(iDim2) == 1) ? (xS) : ((std::abs(iDim2) == 2) ? (yS) : (zS));
        if (iDim3 < 0) zD= (std::abs(iDim3) == 1) ? (nbXS - 1 - xS) : ((std::abs(iDim3) == 2) ? (nbYS - 1 - yS) : (nbZS - 1 - zS));
        if (iDim3 > 0) zD= (std::abs(iDim3) == 1) ? (xS) : ((std::abs(iDim3) == 2) ? (yS) : (zS));
        oRule[0][xD][yD][zD]= iRule[0][xS][yS][zS];
        oRule[1][xD][yD][zD]= iRule[1][xS][yS][zS];
      }
    }
  }

  return oRule;
}


std::array<std::vector<std::vector<std::vector<int>>>, 2> MarkovProcGene::BuildColorSwap(const int iOldColor, const int iNewColor,
                                                                                         const std::array<std::vector<std::vector<std::vector<int>>>, 2>& iRule) {
  if (iRule[0].empty()) throw;
  if (iRule[0][0].empty()) throw;
  if (iRule[0][0][0].empty()) throw;
  int nbXS= (int)iRule[0].size();
  int nbYS= (int)iRule[0][0].size();
  int nbZS= (int)iRule[0][0][0].size();

  std::array<std::vector<std::vector<std::vector<int>>>, 2> oRule({Field::AllocField3D(nbXS, nbYS, nbZS, 0), Field::AllocField3D(nbXS, nbYS, nbZS, 0)});
  for (int xS= 0; xS < nbXS; xS++) {
    for (int yS= 0; yS < nbYS; yS++) {
      for (int zS= 0; zS < nbZS; zS++) {
        oRule[0][xS][yS][zS]= iRule[0][xS][yS][zS];
        oRule[1][xS][yS][zS]= iRule[1][xS][yS][zS];
        if (oRule[0][xS][yS][zS] == iOldColor) oRule[0][xS][yS][zS]= iNewColor;
        if (oRule[1][xS][yS][zS] == iOldColor) oRule[1][xS][yS][zS]= iNewColor;
      }
    }
  }

  return oRule;
}
