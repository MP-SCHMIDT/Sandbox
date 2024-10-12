// Standard lib
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>

// GLUT lib
#include "GL/freeglut.h"

// Trackball lib
#include "Trackball/Camera.hpp"

// Project data
#include "Data.hpp"

// Algo headers
#include "Draw/Colormap.hpp"

// Project classes
#include "AgentSwarmBoid/AgentSwarmBoid.hpp"
#include "AlgoTestEnviro/AlgoTestEnviro.hpp"
#include "BoardGameBotAI/BoardGameBotAI.hpp"
#include "FractalCurvDev/FractalCurvDev.hpp"
#include "FractalElevMap/FractalElevMap.hpp"
#include "ImageExtruMesh/ImageExtruMesh.hpp"
#include "MarkovProcGene/MarkovProcGene.hpp"
#include "MassSpringSyst/MassSpringSyst.hpp"
#include "ParticForceLaw/ParticForceLaw.hpp"
#include "ParticLifeOrga/ParticLifeOrga.hpp"
#include "PosiBasedDynam/PosiBasedDynam.hpp"
#include "SkeletonFolder/SkeletonFolder.hpp"
#include "SpaceTimeWorld/SpaceTimeWorld.hpp"
#include "StringArtOptim/StringArtOptim.hpp"
#include "TerrainErosion/TerrainErosion.hpp"
#include "TestsKernelGPU/TestsKernelGPU.hpp"
#include "WavePropagSimu/WavePropagSimu.hpp"
//#define PRIVATE_RESEARCH_SANDBOX_SUPERSET
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
#include "CompuFluidDyna/CompuFluidDyna.hpp"
#include "NonLinMMABench/NonLinMMABench.hpp"
#include "StructGenOptim/StructGenOptim.hpp"
#endif


// Global variables used by the display
static int windowID;
static int winW, winH;
static int winPosW, winPosH;
static int currentProjectID;
static bool isDarkMode;
static bool isSmoothDraw;
static bool autoScalePlot;
static bool autoScaleScatter;
Camera *cam;

// Global constants used by the display
constexpr int winFPS= 60;             // Target framerate
constexpr int paramPerPage= 40;       // Max number of UI param per page
constexpr int paramLabelNbChar= 16;   // Number of characters per label
constexpr int paramSpaceNbChar= 1;    // Number of character spacing between label and value
constexpr int paramValSignNbChar= 1;  // Number of characters in sign part of the param value
constexpr int paramValInteNbChar= 9;  // Number of characters in integer part of the param value
constexpr int paramValSepaNbChar= 1;  // Number of characters in separator part of the param value
constexpr int paramValFracNbChar= 9;  // Number of characters in fractional part of the param value
constexpr int paramValNbChar= paramValSignNbChar + paramValInteNbChar + paramValSepaNbChar + paramValFracNbChar;
constexpr int plotAreaW= 600;  // Width of the plot area
constexpr int plotAreaH= 120;  // Height of the plot area
constexpr int scatAreaW= 280;  // Width of the scatter area
constexpr int scatAreaH= 280;  // Height of the scatter area
constexpr int winMarginL= 2;   // Margin size on the left of the window
constexpr int winMarginR= 2;   // Margin size on the right of the window
constexpr int winMarginT= 0;   // Margin size on the top of the window
constexpr int winMarginB= 5;   // Margin size on the bottom of the window
constexpr int charHeight= 14;  // Character width in pixels
constexpr int charWidth= 10;   // Character height in pixels
constexpr int textBoxW= 9 * charWidth;
constexpr int textBoxH= charHeight;

// Global variables used by the projects
Data D;
AgentSwarmBoid myAgentSwarmBoid;
AlgoTestEnviro myAlgoTestEnviro;
FractalCurvDev myFractalCurvDev;
FractalElevMap myFractalElevMap;
ImageExtruMesh myImageExtruMesh;
MarkovProcGene myMarkovProcGene;
MassSpringSyst myMassSpringSyst;
ParticForceLaw myParticForceLaw;
ParticLifeOrga myParticLifeOrga;
PosiBasedDynam myPosiBasedDynam;
SkeletonFolder mySkeletonFolder;
SpaceTimeWorld mySpaceTimeWorld;
StringArtOptim myStringArtOptim;
TerrainErosion myTerrainErosion;
TestsKernelGPU myTestsKernelGPU;
BoardGameBotAI myBoardGameBotAI;
WavePropagSimu myWavePropagSimu;
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
CompuFluidDyna myCompuFluidDyna;
NonLinMMABench myNonLinMMABench;
StructGenOptim myStructGenOptim;
#endif

enum ProjectID
{
  AaaaaaaaaaaaaaID,
  AgentSwarmBoidID,
  AlgoTestEnviroID,
  BoardGameBotAIID,
  FractalCurvDevID,
  FractalElevMapID,
  ImageExtruMeshID,
  MarkovProcGeneID,
  MassSpringSystID,
  ParticForceLawID,
  ParticLifeOrgaID,
  PosiBasedDynamID,
  SkeletonFolderID,
  SpaceTimeWorldID,
  StringArtOptimID,
  TerrainErosionID,
  TestsKernelGPUID,
  WavePropagSimuID,
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  CompuFluidDynaID,
  NonLinMMABenchID,
  StructGenOptimID,
#endif
  ZzzzzzzzzzzzzzID,
};

void project_ForceHardInit() {
  D.Plot.clear();
  D.Scatter.clear();
  D.Status.clear();

  D.boxMin= {0.0, 0.0, 0.0};
  D.boxMax= {1.0, 1.0, 1.0};

  int newProjInit= D.UI.empty() ? 1 : 0;
  if (currentProjectID != ProjectID::AgentSwarmBoidID && myAgentSwarmBoid.isActivProj && ++newProjInit) myAgentSwarmBoid.isActivProj= false;
  if (currentProjectID != ProjectID::AlgoTestEnviroID && myAlgoTestEnviro.isActivProj && ++newProjInit) myAlgoTestEnviro.isActivProj= false;
  if (currentProjectID != ProjectID::BoardGameBotAIID && myBoardGameBotAI.isActivProj && ++newProjInit) myBoardGameBotAI.isActivProj= false;
  if (currentProjectID != ProjectID::FractalCurvDevID && myFractalCurvDev.isActivProj && ++newProjInit) myFractalCurvDev.isActivProj= false;
  if (currentProjectID != ProjectID::FractalElevMapID && myFractalElevMap.isActivProj && ++newProjInit) myFractalElevMap.isActivProj= false;
  if (currentProjectID != ProjectID::ImageExtruMeshID && myImageExtruMesh.isActivProj && ++newProjInit) myImageExtruMesh.isActivProj= false;
  if (currentProjectID != ProjectID::MarkovProcGeneID && myMarkovProcGene.isActivProj && ++newProjInit) myMarkovProcGene.isActivProj= false;
  if (currentProjectID != ProjectID::MassSpringSystID && myMassSpringSyst.isActivProj && ++newProjInit) myMassSpringSyst.isActivProj= false;
  if (currentProjectID != ProjectID::ParticForceLawID && myParticForceLaw.isActivProj && ++newProjInit) myParticForceLaw.isActivProj= false;
  if (currentProjectID != ProjectID::ParticLifeOrgaID && myParticLifeOrga.isActivProj && ++newProjInit) myParticLifeOrga.isActivProj= false;
  if (currentProjectID != ProjectID::PosiBasedDynamID && myPosiBasedDynam.isActivProj && ++newProjInit) myPosiBasedDynam.isActivProj= false;
  if (currentProjectID != ProjectID::SkeletonFolderID && mySkeletonFolder.isActivProj && ++newProjInit) mySkeletonFolder.isActivProj= false;
  if (currentProjectID != ProjectID::SpaceTimeWorldID && mySpaceTimeWorld.isActivProj && ++newProjInit) mySpaceTimeWorld.isActivProj= false;
  if (currentProjectID != ProjectID::StringArtOptimID && myStringArtOptim.isActivProj && ++newProjInit) myStringArtOptim.isActivProj= false;
  if (currentProjectID != ProjectID::TerrainErosionID && myTerrainErosion.isActivProj && ++newProjInit) myTerrainErosion.isActivProj= false;
  if (currentProjectID != ProjectID::TestsKernelGPUID && myTestsKernelGPU.isActivProj && ++newProjInit) myTestsKernelGPU.isActivProj= false;
  if (currentProjectID != ProjectID::WavePropagSimuID && myWavePropagSimu.isActivProj && ++newProjInit) myWavePropagSimu.isActivProj= false;
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID != ProjectID::CompuFluidDynaID && myCompuFluidDyna.isActivProj && ++newProjInit) myCompuFluidDyna.isActivProj= false;
  if (currentProjectID != ProjectID::NonLinMMABenchID && myNonLinMMABench.isActivProj && ++newProjInit) myNonLinMMABench.isActivProj= false;
  if (currentProjectID != ProjectID::StructGenOptimID && myStructGenOptim.isActivProj && ++newProjInit) myStructGenOptim.isActivProj= false;
#endif

  if (newProjInit) {
    for (int k= 0; k < 10; k++) {
      D.displayMode[k]= true;
      D.displayModeLabel[k]= "";
    }
  }

  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.SetActiveProject();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.SetActiveProject();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.SetActiveProject();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.SetActiveProject();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.SetActiveProject();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.SetActiveProject();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.SetActiveProject();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.SetActiveProject();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.SetActiveProject();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.SetActiveProject();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.SetActiveProject();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.SetActiveProject();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.SetActiveProject();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.SetActiveProject();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.SetActiveProject();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.SetActiveProject();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.SetActiveProject();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.SetActiveProject();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.SetActiveProject();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.SetActiveProject();
#endif

  if (newProjInit) {
    for (int k= 1; k < 10; k++)
      if (D.displayModeLabel[k].length() == 0)
        D.displayMode[k]= false;
  }

  D.idxFirstParamPageUI= (D.idxFirstParamPageUI < (int)D.UI.size()) ? D.idxFirstParamPageUI : 0;
  D.idxParamUI= (D.idxParamUI < (int)D.UI.size()) ? D.idxParamUI : 0;
}


void project_Refresh() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.Refresh();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.Refresh();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.Refresh();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.Refresh();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.Refresh();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.Refresh();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.Refresh();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.Refresh();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.Refresh();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.Refresh();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.Refresh();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.Refresh();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.Refresh();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.Refresh();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.Refresh();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.Refresh();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.Refresh();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.Refresh();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.Refresh();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.Refresh();
#endif
}


void project_ParamChange() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.ParamChange();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.ParamChange();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.ParamChange();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.ParamChange();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.ParamChange();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.ParamChange();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.ParamChange();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.ParamChange();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.ParamChange();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.ParamChange();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.ParamChange();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.ParamChange();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.ParamChange();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.ParamChange();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.ParamChange();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.ParamChange();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.ParamChange();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.ParamChange();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.ParamChange();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.ParamChange();
#endif
}


void project_KeyPress() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.KeyPress();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.KeyPress();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.KeyPress();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.KeyPress();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.KeyPress();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.KeyPress();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.KeyPress();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.KeyPress();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.KeyPress();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.KeyPress();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.KeyPress();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.KeyPress();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.KeyPress();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.KeyPress();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.KeyPress();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.KeyPress();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.KeyPress();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.KeyPress();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.KeyPress();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.KeyPress();
#endif
}


void project_MousePress() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.MousePress();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.MousePress();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.MousePress();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.MousePress();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.MousePress();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.MousePress();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.MousePress();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.MousePress();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.MousePress();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.MousePress();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.MousePress();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.MousePress();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.MousePress();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.MousePress();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.MousePress();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.MousePress();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.MousePress();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.MousePress();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.MousePress();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.MousePress();
#endif
}

void project_Animate() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.Animate();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.Animate();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.Animate();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.Animate();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.Animate();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.Animate();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.Animate();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.Animate();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.Animate();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.Animate();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.Animate();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.Animate();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.Animate();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.Animate();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.Animate();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.Animate();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.Animate();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.Animate();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.Animate();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.Animate();
#endif
}


void project_Draw() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.Draw();
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.Draw();
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.Draw();
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.Draw();
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.Draw();
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.Draw();
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.Draw();
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.Draw();
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.Draw();
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.Draw();
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.Draw();
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.Draw();
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.Draw();
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.Draw();
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.Draw();
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.Draw();
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.Draw();
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.Draw();
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.Draw();
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.Draw();
#endif
}


void project_QueueSoftRefresh() {
  if (currentProjectID == ProjectID::AgentSwarmBoidID) myAgentSwarmBoid.isRefreshed= false;
  if (currentProjectID == ProjectID::AlgoTestEnviroID) myAlgoTestEnviro.isRefreshed= false;
  if (currentProjectID == ProjectID::BoardGameBotAIID) myBoardGameBotAI.isRefreshed= false;
  if (currentProjectID == ProjectID::FractalCurvDevID) myFractalCurvDev.isRefreshed= false;
  if (currentProjectID == ProjectID::FractalElevMapID) myFractalElevMap.isRefreshed= false;
  if (currentProjectID == ProjectID::ImageExtruMeshID) myImageExtruMesh.isRefreshed= false;
  if (currentProjectID == ProjectID::MarkovProcGeneID) myMarkovProcGene.isRefreshed= false;
  if (currentProjectID == ProjectID::MassSpringSystID) myMassSpringSyst.isRefreshed= false;
  if (currentProjectID == ProjectID::ParticForceLawID) myParticForceLaw.isRefreshed= false;
  if (currentProjectID == ProjectID::ParticLifeOrgaID) myParticLifeOrga.isRefreshed= false;
  if (currentProjectID == ProjectID::PosiBasedDynamID) myPosiBasedDynam.isRefreshed= false;
  if (currentProjectID == ProjectID::SkeletonFolderID) mySkeletonFolder.isRefreshed= false;
  if (currentProjectID == ProjectID::SpaceTimeWorldID) mySpaceTimeWorld.isRefreshed= false;
  if (currentProjectID == ProjectID::StringArtOptimID) myStringArtOptim.isRefreshed= false;
  if (currentProjectID == ProjectID::TerrainErosionID) myTerrainErosion.isRefreshed= false;
  if (currentProjectID == ProjectID::TestsKernelGPUID) myTestsKernelGPU.isRefreshed= false;
  if (currentProjectID == ProjectID::WavePropagSimuID) myWavePropagSimu.isRefreshed= false;
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  if (currentProjectID == ProjectID::CompuFluidDynaID) myCompuFluidDyna.isRefreshed= false;
  if (currentProjectID == ProjectID::NonLinMMABenchID) myNonLinMMABench.isRefreshed= false;
  if (currentProjectID == ProjectID::StructGenOptimID) myStructGenOptim.isRefreshed= false;
#endif
}


// Utility function to save persistent sandbox configuration on disk
void saveConfigWindow() {
  FILE *file= nullptr;
  file= fopen("ConfigWindow.txt", "w");
  if (file != nullptr) {
    fprintf(file, "winPosW winPosH %d %d\n", winPosW, winPosH);
    fprintf(file, "winW winH %d %d\n", winW, winH);
    fclose(file);
  }
}


// Utility function to load persistent sandbox configuration from disk
void loadConfigWindow() {
  FILE *file= nullptr;
  file= fopen("ConfigWindow.txt", "r");
  if (file != nullptr) {
    fscanf(file, "winPosW winPosH %d %d\n", &winPosW, &winPosH);
    fscanf(file, "winW winH %d %d\n", &winW, &winH);
    fclose(file);
  }
}


// Utility function to save persistent project configuration on disk
void saveConfigProject() {
  FILE *file= nullptr;
  file= fopen("ConfigProject.txt", "w");
  if (file != nullptr) {
    fprintf(file, "currentProjectID %d\n", currentProjectID);
    fprintf(file, "nbParam %d\n", (int)D.UI.size());
    for (int idxParam= 0; idxParam < (int)D.UI.size(); idxParam++) {
      fprintf(file, "%s %e\n", D.UI[idxParam].name.c_str(), D.UI[idxParam].D());
    }
    fclose(file);
  }
}


// Utility function to load persistent project configuration from disk
void loadConfigProject() {
  FILE *file= nullptr;
  file= fopen("ConfigProject.txt", "r");
  if (file != nullptr) {
    fscanf(file, "currentProjectID %d\n", &currentProjectID);
    project_ForceHardInit();
    int nbParam= 0;
    fscanf(file, "nbParam %d\n", &nbParam);
    if (nbParam == (int)D.UI.size()) {
      for (int idxParam= 0; idxParam < nbParam; idxParam++) {
        double val= 0.0;
        char name[100];
        if (fscanf(file, "%s %lf\n", name, &val) == 2 && D.UI[idxParam].name == std::string(name)) {
          D.UI[idxParam].Set(val);
        }
        else {
          printf("Invalid parameter count in loaded config file. Using default parameters.\n");
          D.UI.clear();
          project_ForceHardInit();
          break;
        }
      }
    }
    else {
      printf("Invalid parameter names in loaded config file. Using default parameters.\n");
    }
    fclose(file);
  }
}


// Returns the elapsed time since its last call
float elapsed_time() {
  static long last_time= -1;
  float t= 0.0f;
  long current_time= glutGet(GLUT_ELAPSED_TIME);
  if (last_time == -1) {
    last_time= glutGet(GLUT_ELAPSED_TIME);
  }
  else {
    t+= (float)(current_time - last_time) / 1000.0f;
    last_time= current_time;
  }
  return t;
}


// Utility function to draw text
void draw_text(const int w, const int h, const char *text) {
  constexpr float baseHeight= 152.38f;
  constexpr float baseWidth= 104.76f;
  glPushMatrix();
  glTranslatef(float(w), float(h), 0.0f);
  glScalef((float)charWidth / baseWidth, (float)charHeight / baseHeight, 1.0f);
  for (const char *p= text; *p; p++)
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, *p);
  glPopMatrix();
}


void ComputeMouseIn3D(int x, int y) {
  // Set the camera transformation matrix for the scene
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, double(winW) / double(winH), 0.1, 1000.0);

  // Set the world transformation matrix for the scene
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  cam->setWindowSize(float(winW), float(winH));
  glMultMatrixf(cam->getViewMatrix());

  // Unproject the mouse to near and far planes to get 3D position of beginning and end of mouse ray
  double matModelView[16], matProjection[16];
  int viewport[4];
  glGetDoublev(GL_MODELVIEW_MATRIX, matModelView);
  glGetDoublev(GL_PROJECTION_MATRIX, matProjection);
  glGetIntegerv(GL_VIEWPORT, viewport);
  const double winX= (double)x;
  const double winY= viewport[3] - (double)y;
  gluUnProject(winX, winY, 0.0, matModelView, matProjection, viewport, &D.mouseBeg[0], &D.mouseBeg[1], &D.mouseBeg[2]);
  gluUnProject(winX, winY, 1.0, matModelView, matProjection, viewport, &D.mouseEnd[0], &D.mouseEnd[1], &D.mouseEnd[2]);

  // Intersect mouse ray with X Y and Z planes to get 3D projected positions
  const double midX= 0.5 * (D.boxMax[0] + D.boxMin[0]);
  const double midY= 0.5 * (D.boxMax[1] + D.boxMin[1]);
  const double midZ= 0.5 * (D.boxMax[2] + D.boxMin[2]);
  D.mouseProjX[0]= midX;
  D.mouseProjX[1]= D.mouseBeg[1] + (D.mouseEnd[1] - D.mouseBeg[1]) * (midX - D.mouseBeg[0]) / (D.mouseEnd[0] - D.mouseBeg[0]);
  D.mouseProjX[2]= D.mouseBeg[2] + (D.mouseEnd[2] - D.mouseBeg[2]) * (midX - D.mouseBeg[0]) / (D.mouseEnd[0] - D.mouseBeg[0]);
  D.mouseProjY[0]= D.mouseBeg[0] + (D.mouseEnd[0] - D.mouseBeg[0]) * (midY - D.mouseBeg[1]) / (D.mouseEnd[1] - D.mouseBeg[1]);
  D.mouseProjY[1]= midY;
  D.mouseProjY[2]= D.mouseBeg[2] + (D.mouseEnd[2] - D.mouseBeg[2]) * (midY - D.mouseBeg[1]) / (D.mouseEnd[1] - D.mouseBeg[1]);
  D.mouseProjZ[0]= D.mouseBeg[0] + (D.mouseEnd[0] - D.mouseBeg[0]) * (midZ - D.mouseBeg[2]) / (D.mouseEnd[2] - D.mouseBeg[2]);
  D.mouseProjZ[1]= D.mouseBeg[1] + (D.mouseEnd[1] - D.mouseBeg[1]) * (midZ - D.mouseBeg[2]) / (D.mouseEnd[2] - D.mouseBeg[2]);
  D.mouseProjZ[2]= midZ;
}


// Display callback
void callback_display() {
  // Set and clear viewport
  glViewport(0, 0, winW, winH);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set the camera transformation matrix for the scene
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, double(winW) / double(winH), 0.1, 1000.0);

  // Set the world transformation matrix for the scene
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  cam->setWindowSize(float(winW), float(winH));
  const float *viewMat= cam->getViewMatrix();
  glMultMatrixf(viewMat);
  D.camDir= {viewMat[2], viewMat[6], viewMat[10]};

  // Draw the origin unit frame
  if (D.showAxisMode == 1 || D.showAxisMode == 3) {
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();
    glLineWidth(1.0f);

    glPointSize(6.0f);
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();
    glPointSize(1.0f);
  }

  // Draw the reference box
  if (D.showAxisMode == 2 || D.showAxisMode == 3) {
    glColor3f(0.5f, 0.5f, 0.5f);
    glPushMatrix();
    glTranslatef((float)D.boxMin[0], (float)D.boxMin[1], (float)D.boxMin[2]);
    glScalef((float)D.boxMax[0] - (float)D.boxMin[0], (float)D.boxMax[1] - (float)D.boxMin[1], (float)D.boxMax[2] - (float)D.boxMin[2]);
    glTranslatef(0.5f, 0.5f, 0.5f);
    glutWireCube(1.0);
    glPopMatrix();
  }

  // Draw stuff in the scene
  project_Draw();

  // Disable potential clipping before drawing the UI
  glDisable(GL_CLIP_PLANE0);

  // Set the camera transformation matrix for the HUD
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, double(winW), 0.0, double(winH), -1.0, 1.0);

  // Set the world transformation matrix for the HUD
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Add the rectangle backgrounds behind UI elements
  if (isDarkMode) glColor4f(0.4f, 0.4f, 0.4f, 0.5f);
  else glColor4f(0.8f, 0.8f, 0.8f, 0.5f);
  glRecti(0, 0, winMarginL + (paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1.1) * charWidth, winH);
  glRecti(winMarginL + (paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1.1) * charWidth, 0, winW, winMarginB + charHeight);
  glRecti(winMarginL + (paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1.1) * charWidth, winH - winMarginT - plotAreaH - 2.3 * charHeight, winW, winH);

  // Draw the parameter list
  glLineWidth(2.0f);
  for (int k= D.idxFirstParamPageUI; k < std::min((int)D.UI.size(), D.idxFirstParamPageUI + paramPerPage); k++) {
    if (k == D.idxParamUI) {
      if (isDarkMode) glColor3f(0.8f, 0.4f, 0.4f);
      else glColor3f(0.8f, 0.2f, 0.2f);
    }
    else {
      if (isDarkMode) glColor3f(0.8f, 0.8f, 0.8f);
      else glColor3f(0.2f, 0.2f, 0.2f);
    }
    char str[50];
    sprintf(str, "%s %+020.9f", D.UI[k].name.c_str(), D.UI[k].D());  // Format must match paramValNbChar settings
    draw_text(winMarginL, winH - winMarginT - (k - D.idxFirstParamPageUI + 1) * charHeight, str);
    if (k == D.idxParamUI) {
      sprintf(str, "_");
      draw_text(winMarginL + (paramLabelNbChar + paramSpaceNbChar + D.idxCursorUI) * charWidth,
                winH - winMarginT - (k - D.idxFirstParamPageUI + 1) * charHeight, str);
      draw_text(winMarginL + (paramLabelNbChar + paramSpaceNbChar + D.idxCursorUI) * charWidth,
                winH - winMarginT - (k - D.idxFirstParamPageUI) * charHeight, str);
    }
  }
  glLineWidth(1.0f);

  // Draw the parameter list scroll bar
  if (!D.UI.empty()) {
    if (isDarkMode) glColor3f(0.8f, 0.8f, 0.8f);
    else glColor3f(0.2f, 0.2f, 0.2f);
    glBegin(GL_LINES);
    glVertex3i((paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1) * charWidth,
               winH - winMarginT, 0);
    glVertex3i((paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1) * charWidth,
               winH - winMarginT - textBoxH * paramPerPage, 0);
    glEnd();
    glLineWidth(4.0f);
    glBegin(GL_LINES);
    const float barBeg= (float)D.idxFirstParamPageUI / (float)D.UI.size();
    const float barEnd= (float)std::min(D.idxFirstParamPageUI + paramPerPage, (int)D.UI.size()) / (float)D.UI.size();
    glVertex3i((paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1) * charWidth,
               winH - winMarginT - textBoxH * paramPerPage * barBeg, 0);
    glVertex3i((paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 1) * charWidth,
               winH - winMarginT - textBoxH * paramPerPage * barEnd, 0);
    glEnd();
    glLineWidth(1.0f);
  }

  // Draw the 2D plot
  if (!D.Plot.empty()) {
    double valMin= 0.0, valMax= 1.0;
    for (int k0= 0; k0 < (int)D.Plot.size(); k0++) {
      if (D.Plot[k0].val.empty()) continue;

      // Set the color
      float r, g, b;
      Colormap::RatioToRainbow(float(k0) / (float)std::max((int)D.Plot.size() - 1, 1), r, g, b);
      glColor3f(r, g, b);

      // Find the min max range for vertical scaling
      if (autoScalePlot && (k0 == 0 || !D.Plot[k0].isSameRange)) {
        valMin= std::numeric_limits<double>::max();
        valMax= std::numeric_limits<double>::lowest();
        // Also check following plots if using the same range
        for (int k1= k0; k1 < (int)D.Plot.size(); k1++) {
          if (k1 > k0 && !D.Plot[k1].isSameRange) break;
          for (double valCur : D.Plot[k1].val) {
            if (valMin > valCur) valMin= valCur;
            if (valMax < valCur) valMax= valCur;
          }
          if (D.Plot[k1].isSymmetric) {
            valMax= std::max(std::abs(valMin), std::abs(valMax));
            valMin= -valMax;
          }
        }
      }

      // Draw the text for legend, min max start end values
      glLineWidth(2.0f);
      char str[50];
      sprintf(str, "N:%d", int(D.Plot[k0].val.size()));  // Count
      draw_text(winW - winMarginR - plotAreaW - 5 * textBoxW,
                winH - winMarginT - textBoxH * (k0 + 2), str);
      strcpy(str, D.Plot[k0].name.c_str());  // Name
      draw_text(winW - winMarginR - plotAreaW - 4 * textBoxW,
                winH - winMarginT - textBoxH * (k0 + 2), str);
      if (D.Plot[k0].isLog) {  // Log flag
        sprintf(str, "log");
        draw_text(winW - winMarginR - plotAreaW - 2.5 * textBoxW,
                  winH - winMarginT - textBoxH * (k0 + 2), str);
      }
      const int spacingMinMaxVal= (plotAreaW - textBoxW) / std::max(int(D.Plot.size() - 1), 1);
      sprintf(str, "%+.2e", valMax);  // Max value
      draw_text(winW - winMarginR - textBoxW - plotAreaW + k0 * spacingMinMaxVal,
                winH - winMarginT - textBoxH, str);
      sprintf(str, "%+.2e", valMin);  // Min value
      draw_text(winW - winMarginR - textBoxW - plotAreaW + k0 * spacingMinMaxVal,
                winH - winMarginT - plotAreaH - 2 * textBoxH, str);
      sprintf(str, "%+.2e", D.Plot[k0].val[0]);  // Start value
      draw_text(winW - winMarginR - plotAreaW - 2 * textBoxW,
                winH - winMarginT - textBoxH - textBoxH * k0 - textBoxH, str);
      sprintf(str, "%+.2e", D.Plot[k0].val[D.Plot[k0].val.size() - 1]);  // End value
      draw_text(winW - winMarginR - textBoxW,
                winH - winMarginT - textBoxH - textBoxH * k0 - textBoxH, str);
      glLineWidth(1.0f);

      // Draw the zero axes
      if (valMax > valMin && !D.Plot[k0].isLog) {
        double valScaled= -valMin / (valMax - valMin);
        if (valScaled >= 0.0 && valScaled <= 1.0) {
          glBegin(GL_LINES);
          glVertex3i(winW - winMarginR - plotAreaW - textBoxW,
                     winH - winMarginT - plotAreaH - textBoxH + plotAreaH * valScaled, 0);
          glVertex3i(winW - winMarginR - plotAreaW - textBoxW + plotAreaW,
                     winH - winMarginT - plotAreaH - textBoxH + plotAreaH * valScaled, 0);
          glEnd();
        }
      }

      // Draw the plot curves and markers
      glLineWidth(2.0f);
      glPointSize(3.0f);
      for (int mode= 0; mode < (D.Plot[k0].showPoints ? 2 : 1); mode++) {
        if (mode == 0) glBegin(GL_LINE_STRIP);
        if (mode == 1) glBegin(GL_POINTS);
        for (int k1= 0; k1 < (int)D.Plot[k0].val.size(); k1++) {
          double valScaled;
          if (valMax - valMin == 0.0) valScaled= 0.0;
          else if (!D.Plot[k0].isLog) valScaled= (D.Plot[k0].val[k1] - valMin) / (valMax - valMin);
          else if (D.Plot[k0].val[k1] <= 0.0) valScaled= 0.0;
          else if (valMin < 0.0) valScaled= 1.0;
          else valScaled= (std::log10(D.Plot[k0].val[k1]) - std::log10(valMin)) / (std::log10(valMax) - std::log10(valMin));
          glVertex3i(winW - winMarginR - plotAreaW - textBoxW + plotAreaW * k1 / std::max((int)D.Plot[k0].val.size() - 1, 1),
                     winH - winMarginT - plotAreaH - textBoxH + plotAreaH * valScaled, 0);
        }
        glEnd();
      }
      glLineWidth(1.0f);
      glPointSize(1.0f);
    }
  }

  // Draw the 2D scatter
  if (!D.Scatter.empty()) {
    // Find the min max range for scaling
    double valMinX= -1.0;
    double valMinY= -1.0;
    double valMaxX= 1.0;
    double valMaxY= 1.0;
    if (autoScaleScatter) {
      valMinX= std::numeric_limits<double>::max();
      valMinY= std::numeric_limits<double>::max();
      valMaxX= std::numeric_limits<double>::lowest();
      valMaxY= std::numeric_limits<double>::lowest();
      for (int k0= 0; k0 < int(D.Scatter.size()); k0++) {
        for (int k1= 0; k1 < int(D.Scatter[k0].val.size()); k1++) {
          if (valMinX > D.Scatter[k0].val[k1][0]) valMinX= D.Scatter[k0].val[k1][0];
          if (valMinY > D.Scatter[k0].val[k1][1]) valMinY= D.Scatter[k0].val[k1][1];
          if (valMaxX < D.Scatter[k0].val[k1][0]) valMaxX= D.Scatter[k0].val[k1][0];
          if (valMaxY < D.Scatter[k0].val[k1][1]) valMaxY= D.Scatter[k0].val[k1][1];
        }
      }
    }

    // Draw the axes and zero lines
    if (isDarkMode) glColor3f(0.8f, 0.8f, 0.8f);
    else glColor3f(0.2f, 0.2f, 0.2f);
    glBegin(GL_LINE_LOOP);
    glVertex3i(winMarginL + textBoxW, winMarginB + scatAreaH + 3 * textBoxH, 0);
    glVertex3i(winMarginL + textBoxW, winMarginB + 3 * textBoxH, 0);
    glVertex3i(winMarginL + textBoxW + scatAreaW, winMarginB + 3 * textBoxH, 0);
    glVertex3i(winMarginL + textBoxW + scatAreaW, winMarginB + scatAreaH + 3 * textBoxH, 0);
    glEnd();
    if (valMinX < 0.0 && valMaxX > 0.0) {
      int offsetW= ((0.0 - valMinX) / (valMaxX - valMinX)) * scatAreaW;
      glBegin(GL_LINES);
      glVertex3i(winMarginL + textBoxW + offsetW, winMarginB + scatAreaH + 3 * textBoxH, 0);
      glVertex3i(winMarginL + textBoxW + offsetW, winMarginB + 3 * textBoxH, 0);
      glEnd();
    }
    if (valMinY < 0.0 && valMaxY > 0.0) {
      int offsetH= ((0.0 - valMinY) / (valMaxY - valMinY)) * scatAreaH;
      glBegin(GL_LINES);
      glVertex3i(winMarginL + textBoxW, winMarginB + offsetH + 3 * textBoxH, 0);
      glVertex3i(winMarginL + textBoxW + scatAreaW, winMarginB + offsetH + 3 * textBoxH, 0);
      glEnd();
    }

    // Draw the min max values
    char str[50];
    glLineWidth(2.0f);
    sprintf(str, "%+.2e", valMinX);
    draw_text(winMarginL + textBoxW, winMarginB + 2 * textBoxH, str);
    sprintf(str, "%+.2e", valMaxX);
    draw_text(winMarginL + textBoxW + scatAreaW - textBoxW, winMarginB + 2 * textBoxH, str);
    sprintf(str, "%+.2e", valMinY);
    draw_text(winMarginL, winMarginB + 3 * textBoxH, str);
    sprintf(str, "%+.2e", valMaxY);
    draw_text(winMarginL, winMarginB + 3 * textBoxH + scatAreaH - textBoxH, str);
    glLineWidth(1.0f);

    // Draw the scatter legend and points
    glLineWidth(2.0f);
    glPointSize(3.0f);
    for (int k0= 0; k0 < (int)D.Scatter.size(); k0++) {
      // Set the color
      float r, g, b;
      Colormap::RatioToRainbow(float(k0) / (float)std::max((int)D.Scatter.size() - 1, 1), r, g, b);
      glColor3f(r, g, b);

      // Draw the text for legend and point count
      strcpy(str, D.Scatter[k0].name.c_str());
      draw_text(0, scatAreaH - 2 * k0 * textBoxH, str);
      sprintf(str, "N:%d", int(D.Scatter[k0].val.size()));
      draw_text(0, scatAreaH - (2 * k0 + 1) * textBoxH, str);

      // Draw the points
      glBegin(GL_POINTS);
      for (int k1= 0; k1 < int(D.Scatter[k0].val.size()); k1++) {
        const double relPosX= (D.Scatter[k0].val[k1][0] - valMinX) / (valMaxX - valMinX);
        const double relPosY= (D.Scatter[k0].val[k1][1] - valMinY) / (valMaxY - valMinY);
        glVertex3i(winMarginL + textBoxW + (int)std::round((double)scatAreaW * relPosX),
                   winMarginB + 3 * textBoxH + (int)std::round((double)scatAreaH * relPosY), 0);
      }
      glEnd();
    }
    glLineWidth(1.0f);
    glPointSize(1.0f);
  }

  // Draw the frame time
  {
    glLineWidth(2.0f);
    if (D.playAnimation) {
      if (isDarkMode) glColor3f(0.8f, 0.4f, 0.4f);
      else glColor3f(0.8f, 0.2f, 0.2f);
    }
    else {
      if (isDarkMode) glColor3f(0.8f, 0.8f, 0.8f);
      else glColor3f(0.2f, 0.2f, 0.2f);
    }
    char str[50];
    sprintf(str, "%.3fs", elapsed_time());
    draw_text(winMarginL, winMarginB, str);
    glLineWidth(1.0f);
  }

  // Draw the display mode labels
  if (currentProjectID != 0) {
    glLineWidth(2.0f);
    for (int k= 1; k < 10; k++) {
      if (D.displayMode[k]) {
        if (isDarkMode) glColor3f(0.8f, 0.8f, 0.8f);
        else glColor3f(0.2f, 0.2f, 0.2f);
      }
      else {
        glColor3f(0.5f, 0.5f, 0.5f);
      }
      char str[50];
      sprintf(str, "%d:%s", k, D.displayModeLabel[k].c_str());
      draw_text((paramLabelNbChar + paramSpaceNbChar + paramValNbChar + 2) * charWidth, winH - winMarginT - k * charHeight, str);
    }
    glLineWidth(1.0f);
  }

  // Draw the status text
  {
    glLineWidth(2.0f);
    int cumulLen= 10;
    for (std::string text : D.Status) {
      if (isDarkMode) glColor3f(0.8f, 0.8f, 0.8f);
      else glColor3f(0.2f, 0.2f, 0.2f);
      char str[50];
      sprintf(str, "%s", text.c_str());
      draw_text(winMarginL + cumulLen * charWidth, winMarginB, str);
      cumulLen+= text.length() + 3;
    }
    glLineWidth(1.0f);
  }

  // Commit the draw
  glutSwapBuffers();
}


// Timer program interruption callback
void callback_timer(int v) {
  // Compute animations
  if (D.playAnimation || D.stepAnimation) {
    project_Animate();
    glutPostRedisplay();
    D.stepAnimation= false;
  }

  glutTimerFunc(1000 / winFPS, callback_timer, v);
}


// Window reshape callback
void callback_reshape(int width, int height) {
  winW= width;
  winH= height;
  glutPostRedisplay();
}


// Keyboard interruption callback
void callback_keyboard(unsigned char key, int x, int y) {
  // // Display pressed key code
  // printf("%c  %d\n", key, key);
  D.keyIsShift= (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
  D.keyIsCtrl= (glutGetModifiers() & GLUT_ACTIVE_CTRL);
  D.keyIsAlt= (glutGetModifiers() & GLUT_ACTIVE_ALT);

  if (key == 27) {
    glutDestroyWindow(windowID);
    exit(EXIT_SUCCESS);
  }
  else if (key == '\\') glutPositionWindow(0, 0);
  else if (key == ' ') D.playAnimation= !D.playAnimation;
  else if (key == '.') D.stepAnimation= !D.stepAnimation;
  else if (key == '\b') {
    if (!std::isnan(D.UI[D.idxParamUI].D())) {
      D.UI[D.idxParamUI].Set(0.0);
      project_ParamChange();
    }
  }
  else if (key == '\t' && !D.keyIsShift) D.idxFirstParamPageUI= (D.idxFirstParamPageUI + paramPerPage < (int)D.UI.size()) ? (D.idxFirstParamPageUI + paramPerPage) : (0);
  else if (key == '\t' && D.keyIsShift) D.idxFirstParamPageUI= (D.idxFirstParamPageUI - paramPerPage >= 0) ? (D.idxFirstParamPageUI - paramPerPage) : (((int)D.UI.size() / paramPerPage) * paramPerPage);
  else if (key >= '1' && key <= '9') D.displayMode[key - '0']= !D.displayMode[key - '0'];
  else if (key == '0') D.showAxisMode= (D.showAxisMode + 1) % 4;
  else if (key == '=') {
    D.Plot.clear();
    D.Scatter.clear();
    D.Status.clear();
  }
  else if (key == ',') project_ForceHardInit();
  else if (key == '/') project_QueueSoftRefresh();
  else if ((key >= 'A' && key <= 'Z') || (key >= 'a' && key <= 'z')) {
    ComputeMouseIn3D(x, y);
    D.keyLetterUpperCase= (key >= 'a' && key <= 'z') ? (key - ('a' - 'A')) : (key);
    project_KeyPress();
    D.keyLetterUpperCase= 0;
  }

  // Compute refresh
  if (key != ' ' && key != '.') {
    project_Refresh();
    glutPostRedisplay();
  }
}


// Keyboard special interruption callback
void callback_keyboard_special(int key, int x, int y) {
  (void)x;  // Disable warning unused variable
  (void)y;  // Disable warning unused variable
  D.keyIsShift= (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
  D.keyIsCtrl= (glutGetModifiers() & GLUT_ACTIVE_CTRL);
  D.keyIsAlt= (glutGetModifiers() & GLUT_ACTIVE_ALT);

  if (D.UI.empty()) return;

  // Handle arrow key to switch selected parameter
  if (key == GLUT_KEY_UP || key == GLUT_KEY_DOWN) {
    if (D.keyIsShift) {
      if (key == GLUT_KEY_UP) D.idxParamUI= (D.idxParamUI - 5 + int(D.UI.size())) % int(D.UI.size());
      if (key == GLUT_KEY_DOWN) D.idxParamUI= (D.idxParamUI + 5) % int(D.UI.size());
      D.idxFirstParamPageUI= paramPerPage * (D.idxParamUI / paramPerPage);
    }
    else {
      if (key == GLUT_KEY_UP) D.idxParamUI= (D.idxParamUI - 1 + int(D.UI.size())) % int(D.UI.size());
      if (key == GLUT_KEY_DOWN) D.idxParamUI= (D.idxParamUI + 1) % int(D.UI.size());
      D.idxFirstParamPageUI= paramPerPage * (D.idxParamUI / paramPerPage);
    }
    glutPostRedisplay();
  }

  // Handle arrow key to change the selected parameter value
  if (key == GLUT_KEY_LEFT || key == GLUT_KEY_RIGHT) {
    if (D.keyIsShift) {
      if (key == GLUT_KEY_LEFT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() / 2.0);
      if (key == GLUT_KEY_RIGHT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() * 2.0);
    }
    else if (D.keyIsCtrl) {
      if (key == GLUT_KEY_LEFT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() / (1.0 + 1.0 / 16.0));
      if (key == GLUT_KEY_RIGHT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() * (1.0 + 1.0 / 16.0));
    }
    else if (D.keyIsAlt) {
      if (key == GLUT_KEY_LEFT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() - 1.0 / 16.0);
      if (key == GLUT_KEY_RIGHT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() + 1.0 / 16.0);
    }
    else {
      if (key == GLUT_KEY_LEFT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() - 1.0);
      if (key == GLUT_KEY_RIGHT) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() + 1.0);
    }
    // Compute refresh
    project_ParamChange();
    project_Refresh();
    glutPostRedisplay();
  }
}


// Mouse click interruption callback
void callback_mouse_click(int button, int state, int x, int y) {
  D.keyIsShift= (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
  D.keyIsCtrl= (glutGetModifiers() & GLUT_ACTIVE_CTRL);
  D.keyIsAlt= (glutGetModifiers() & GLUT_ACTIVE_ALT);

  cam->setCurrentMousePos(float(x), float(y));

  if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON && !D.keyIsShift && !D.keyIsCtrl) cam->beginRotate();
  if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON && !D.keyIsShift && D.keyIsCtrl) cam->beginZoom();
  if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON && D.keyIsShift && !D.keyIsCtrl) cam->beginPan();
  if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) cam->endRotate();
  if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) cam->endZoom();
  if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) cam->endPan();

  if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON) {
    ComputeMouseIn3D(x, y);
    D.mouseMiddleButtonState= 1;
    project_MousePress();
    D.mouseMiddleButtonState= 2;
  }
  if (state == GLUT_UP && button == GLUT_MIDDLE_BUTTON) {
    ComputeMouseIn3D(x, y);
    D.mouseMiddleButtonState= 3;
    project_MousePress();
    D.mouseMiddleButtonState= 0;
  }

  if (state == GLUT_UP && (button == 3 || button == 4)) {
    if (!D.UI.empty()) {
      if (x > (paramLabelNbChar + paramSpaceNbChar) * charWidth) {
        if (x < (paramLabelNbChar + paramSpaceNbChar + paramValNbChar) * charWidth) {
          // Change the targeted symbol if any
          bool hasParamChanged= true;
          if (D.idxCursorUI < paramValSignNbChar) {
            D.UI[D.idxParamUI].Set(-D.UI[D.idxParamUI].D());
          }
          else if (D.idxCursorUI >= paramValSignNbChar && D.idxCursorUI < paramValSignNbChar + paramValInteNbChar) {
            if (button == 3) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() + std::pow(10.0, double(paramValInteNbChar - D.idxCursorUI)));
            if (button == 4) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() - std::pow(10.0, double(paramValInteNbChar - D.idxCursorUI)));
          }
          else if (D.idxCursorUI >= paramValSignNbChar + paramValInteNbChar + paramValSepaNbChar && D.idxCursorUI < paramValNbChar) {
            if (button == 3) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() + std::pow(10.0, double(paramValInteNbChar + paramValSepaNbChar - D.idxCursorUI)));
            if (button == 4) D.UI[D.idxParamUI].Set(D.UI[D.idxParamUI].D() - std::pow(10.0, double(paramValInteNbChar + paramValSepaNbChar - D.idxCursorUI)));
          }
          else {
            hasParamChanged= false;
          }
          // Compute refresh
          if (hasParamChanged) {
            project_ParamChange();
            project_Refresh();
          }
        }
      }
    }
  }

  glutPostRedisplay();
}


// Mouse active motion interruption callback
void callback_mouse_motion(int x, int y) {
  cam->setCurrentMousePos(float(x), float(y));

  if (D.mouseMiddleButtonState == 2) {
    D.keyIsShift= (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
    D.keyIsCtrl= (glutGetModifiers() & GLUT_ACTIVE_CTRL);
    D.keyIsAlt= (glutGetModifiers() & GLUT_ACTIVE_ALT);
    ComputeMouseIn3D(x, y);
    project_MousePress();
  }

  glutPostRedisplay();
}


// Mouse passive motion interruption callback
void callback_passive_mouse_motion(int x, int y) {
  const int prevParamIdx= D.idxParamUI;
  const int prevCursorIdx= D.idxCursorUI;
  if (x > winMarginL + (paramLabelNbChar + paramSpaceNbChar) * charWidth &&
      x < winMarginL + (paramLabelNbChar + paramSpaceNbChar + paramValNbChar) * charWidth) {
    if (y > winMarginT &&
        y < winMarginT + std::min((int)D.UI.size() - D.idxFirstParamPageUI, paramPerPage) * charHeight) {
      D.idxParamUI= D.idxFirstParamPageUI + (y - winMarginT) / charHeight;
      D.idxCursorUI= std::min(std::max((x - winMarginL - (paramLabelNbChar + paramSpaceNbChar) * charWidth) / charWidth, 0), paramValNbChar);
    }
  }

  if (D.idxParamUI != prevParamIdx || D.idxCursorUI != prevCursorIdx)
    glutPostRedisplay();
}


// Menu selection callback
void callback_menu(int num) {
  // Reset or activate the selected project
  if (num > ProjectID::AaaaaaaaaaaaaaID && num < ProjectID::ZzzzzzzzzzzzzzID) {
    currentProjectID= num;
    D.UI.clear();
    project_ForceHardInit();
  }
  if (num <= -1 && num >= -6) {
    delete cam;
    cam= new Camera;
    std::array<double, 3> center= {0.5 * (D.boxMin[0] + D.boxMax[0]), 0.5 * (D.boxMin[1] + D.boxMax[1]), 0.5 * (D.boxMin[2] + D.boxMax[2])};
    double diag= std::sqrt(std::pow((D.boxMax[0] - D.boxMin[0]), 2.0) + std::pow((D.boxMax[1] - D.boxMin[1]), 2.0) + std::pow((D.boxMax[2] - D.boxMin[2]), 2.0));
    if (num == -1) cam->setEye(center[0] + diag, center[1], center[2]);
    if (num == -2) cam->setEye(center[0] - diag, center[1], center[2]);
    if (num == -3) cam->setEye(center[0], center[1] + diag, center[2]);
    if (num == -4) cam->setEye(center[0], center[1] - diag, center[2]);
    if (num == -5) cam->setEye(center[0] + 1.e-6 * diag, center[1], center[2] + diag);
    if (num == -6) cam->setEye(center[0] + 1.e-6 * diag, center[1], center[2] - diag);
    cam->setCenter(center[0], center[1], center[2]);
  }
  // Toggle dark mode display
  if (num == -10) {
    isDarkMode= !isDarkMode;
    if (isDarkMode) glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    else glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  }
  // Toggle smooth drawing option
  if (num == -11) {
    isSmoothDraw= !isSmoothDraw;
    if (isSmoothDraw) {
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
    }
    else {
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_LINE_SMOOTH);
    }
  }
  // Toggle auto scale plot
  if (num == -12) {
    autoScalePlot= !autoScalePlot;
  }
  // Toggle auto scale scatter
  if (num == -13) {
    autoScaleScatter= !autoScaleScatter;
  }
  // Save sandbox settings
  if (num == -20) {
    winPosW= glutGet((GLenum)GLUT_WINDOW_X);
    winPosH= glutGet((GLenum)GLUT_WINDOW_Y);
    saveConfigWindow();
  }
  // Save project parameters
  if (num == -21) {
    saveConfigProject();
  }
  // Compute refresh
  project_Refresh();

  glutPostRedisplay();
}


// Menu initialization
void init_menu() {
  // Create menu tree starting from leaves
  const int menuDisplay= glutCreateMenu(callback_menu);
  glutAddMenuEntry("View X+", -1);
  glutAddMenuEntry("View X-", -2);
  glutAddMenuEntry("View Y+", -3);
  glutAddMenuEntry("View Y-", -4);
  glutAddMenuEntry("View Z+", -5);
  glutAddMenuEntry("View Z-", -6);
  glutAddMenuEntry("Dark mode", -10);
  glutAddMenuEntry("Smooth draw", -11);
  glutAddMenuEntry("Auto scale plot", -12);
  glutAddMenuEntry("Auto scale scatter", -13);
  const int menuProject= glutCreateMenu(callback_menu);
  glutAddMenuEntry("AgentSwarmBoid", ProjectID::AgentSwarmBoidID);
  glutAddMenuEntry("AlgoTestEnviro", ProjectID::AlgoTestEnviroID);
  glutAddMenuEntry("BoardGameBotAI", ProjectID::BoardGameBotAIID);
  glutAddMenuEntry("FractalCurvDev", ProjectID::FractalCurvDevID);
  glutAddMenuEntry("FractalElevMap", ProjectID::FractalElevMapID);
  glutAddMenuEntry("ImageExtruMesh", ProjectID::ImageExtruMeshID);
  glutAddMenuEntry("MarkovProcGene", ProjectID::MarkovProcGeneID);
  glutAddMenuEntry("MassSpringSyst", ProjectID::MassSpringSystID);
  glutAddMenuEntry("ParticForceLaw", ProjectID::ParticForceLawID);
  glutAddMenuEntry("ParticLifeOrga", ProjectID::ParticLifeOrgaID);
  glutAddMenuEntry("PosiBasedDynam", ProjectID::PosiBasedDynamID);
  glutAddMenuEntry("SkeletonFolder", ProjectID::SkeletonFolderID);
  glutAddMenuEntry("SpaceTimeWorld", ProjectID::SpaceTimeWorldID);
  glutAddMenuEntry("StringArtOptim", ProjectID::StringArtOptimID);
  glutAddMenuEntry("TerrainErosion", ProjectID::TerrainErosionID);
  glutAddMenuEntry("TestsKernelGPU", ProjectID::TestsKernelGPUID);
  glutAddMenuEntry("WavePropagSimu", ProjectID::WavePropagSimuID);
#ifdef PRIVATE_RESEARCH_SANDBOX_SUPERSET
  glutAddMenuEntry("CompuFluidDyna", ProjectID::CompuFluidDynaID);
  glutAddMenuEntry("NonLinMMABench", ProjectID::NonLinMMABenchID);
  glutAddMenuEntry("StructGenOptim", ProjectID::StructGenOptimID);
#endif
  const int menuSave= glutCreateMenu(callback_menu);
  glutAddMenuEntry("Save settings", -20);
  glutAddMenuEntry("Save parameters", -21);
  glutCreateMenu(callback_menu);
  glutAddSubMenu("Display", menuDisplay);
  glutAddSubMenu("Project", menuProject);
  glutAddSubMenu("Save", menuSave);

  // Attach menu to click
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}


// OpenGL initialization
void init_GL() {
  isDarkMode= true;
  isSmoothDraw= false;
  autoScalePlot= true;
  autoScaleScatter= true;

  // Set background color
  glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

  // Define light properties
  constexpr float light0_ambient[4]= {0.1f, 0.1f, 0.1f, 1.0f};
  constexpr float light0_diffuse[4]= {0.7f, 0.7f, 0.7f, 1.0f};
  constexpr float light0_specular[4]= {0.8f, 0.8f, 0.8f, 1.0f};
  constexpr float position0[4]= {4.0f, 4.0f, 4.0f, 1.0f};
  constexpr float global_ambient[4]= {0.2f, 0.2f, 0.2f, 1.0f};
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, position0);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

  // Define material shader properties
  constexpr float mat_ambiant[4]= {0.2f, 0.2f, 0.2f, 1.0f};
  constexpr float mat_diffuse[4]= {0.8f, 0.8f, 0.8f, 1.0f};
  constexpr float mat_specular[4]= {0.2f, 0.2f, 0.2f, 1.0f};
  constexpr float mat_emission[4]= {0.0f, 0.0f, 0.0f, 1.0f};
  constexpr float mat_shininess[1]= {64.0f};
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambiant);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

  // Misc
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glShadeModel(GL_SMOOTH);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHT0);
}


// Scene initialization
void init_scene() {
  // Initialize pseudo random number generator
  srand(time(0));

  // Initialize camera and arcball positions
  cam= new Camera;
  cam->setEye(2.5f, 0.3f, 0.6f);
  cam->setCenter(0.5f, 0.3f, 0.6f);
}


// Main function
int main(int argc, char *argv[]) {
  // Load window settings or use default values
  winW= 1600;
  winH= 900;
  winPosW= -1;
  winPosH= -1;
  loadConfigWindow();

  // Window creation
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(winW, winH);
  glutInitWindowPosition(winPosW, winPosH);
  windowID= glutCreateWindow("Sandbox");

  // World initialization
  init_menu();
  init_GL();
  init_scene();

  // Callbacks affectation
  glutDisplayFunc(&callback_display);
  glutReshapeFunc(&callback_reshape);
  glutKeyboardFunc(&callback_keyboard);
  glutSpecialFunc(&callback_keyboard_special);
  glutMouseFunc(&callback_mouse_click);
  glutMotionFunc(&callback_mouse_motion);
  glutPassiveMotionFunc(&callback_passive_mouse_motion);
  glutTimerFunc(100, &callback_timer, 0);

  // Compute refresh
  currentProjectID= 0;
  loadConfigProject();
  project_Refresh();

  // Start refresh loop
  glutMainLoop();

  return 0;
}
