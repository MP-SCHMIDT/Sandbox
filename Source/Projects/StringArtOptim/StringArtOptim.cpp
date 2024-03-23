#include "StringArtOptim.hpp"


// Standard lib
#include <cmath>
#include <numbers>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "Draw/Colormap.hpp"
#include "FileIO/FileInput.hpp"
#include "Geom/Bresenham.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


// Constructor
StringArtOptim::StringArtOptim() {
  isActivProj= false;
  isAllocated= false;
  isRefreshed= false;
}


// Initialize Project UI parameters
void StringArtOptim::SetActiveProject() {
  if (!isActivProj) {
    D.UI.clear();
    D.UI.push_back(ParamUI("ImageID_________", 2));
    D.UI.push_back(ParamUI("ImageSizeW______", 256));
    D.UI.push_back(ParamUI("ImageSizeH______", 256));
    D.UI.push_back(ParamUI("PegLayout_______", 1));
    D.UI.push_back(ParamUI("PegNumber_______", 256));
    D.UI.push_back(ParamUI("ColorsAdd_______", 2));
    D.UI.push_back(ParamUI("ColorsSub_______", 0));
    D.UI.push_back(ParamUI("ColorsNormalize_", 0));
    D.UI.push_back(ParamUI("StepCount_______", 1));
    D.UI.push_back(ParamUI("SingleLine______", -0.5));
    D.UI.push_back(ParamUI("BlendMode_______", -0.5));
    D.UI.push_back(ParamUI("CoeffColor______", 0.1));
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
bool StringArtOptim::CheckAlloc() {
  if (D.UI[ImageID_________].hasChanged()) isAllocated= false;
  if (D.UI[ImageSizeW______].hasChanged()) isAllocated= false;
  if (D.UI[ImageSizeH______].hasChanged()) isAllocated= false;
  if (D.UI[PegLayout_______].hasChanged()) isAllocated= false;
  if (D.UI[PegNumber_______].hasChanged()) isAllocated= false;
  if (D.UI[ColorsAdd_______].hasChanged()) isAllocated= false;
  if (D.UI[ColorsSub_______].hasChanged()) isAllocated= false;
  if (D.UI[ColorsNormalize_].hasChanged()) isAllocated= false;
  return isAllocated;
}


// Check if parameter changes should trigger a refresh
bool StringArtOptim::CheckRefresh() {
  return isRefreshed;
}


// Allocate the project data
void StringArtOptim::Allocate() {
  if (!isActivProj) return;
  if (CheckAlloc()) return;
  isRefreshed= false;
  isAllocated= true;

  // Get UI parameters
  nW= std::max(D.UI[ImageSizeW______].I(), 1);
  nH= std::max(D.UI[ImageSizeH______].I(), 1);
}


// Refresh the project
void StringArtOptim::Refresh() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (CheckRefresh()) return;
  isRefreshed= true;

  // Reset plot
  D.Plot.clear();
  D.Scatter.clear();

  // Initialize images
  ImRef= Field::AllocField2D(nW, nH, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  ImCur= Field::AllocField2D(nW, nH, Vec::Vec3<float>(0.0f, 0.0f, 0.0f));
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (D.UI[ImageID_________].I() == 0) FileInput::LoadImageBMPFile("./FileInput/Images/SimpleSmile.bmp", imageRGBA, false);
  if (D.UI[ImageID_________].I() == 1) FileInput::LoadImageBMPFile("./FileInput/Images/MonaLisa.bmp", imageRGBA, false);
  if (D.UI[ImageID_________].I() == 2) FileInput::LoadImageBMPFile("./FileInput/Images/PillarsCreation.bmp", imageRGBA, false);
  if (D.UI[ImageID_________].I() == 3) FileInput::LoadImageBMPFile("./FileInput/Images/Butterfly.bmp", imageRGBA, false);
  if (D.UI[ImageID_________].I() == 4) FileInput::LoadImageBMPFile("./FileInput/Images/AlbertArt.bmp", imageRGBA, false);
  if (D.UI[ImageID_________].I() == 5) FileInput::LoadImageBMPFile("./FileInput/Images/DeepField.bmp", imageRGBA, false);
  if (D.UI[ImageID_________].I() == 6) FileInput::LoadImageBMPFile("./FileInput/Images/Eye.bmp", imageRGBA, false);
  for (int w= 0; w < nW; w++) {
    for (int h= 0; h < nH; h++) {
      if (!imageRGBA.empty()) {
        const float posW= (float)(imageRGBA.size() - 1) * ((float)w + 0.5f) / (float)nW;
        const float posH= (float)(imageRGBA[0].size() - 1) * ((float)h + 0.5f) / (float)nH;
        const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
        const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
        const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
        ImRef[w][h].set(colRGBA[0], colRGBA[1], colRGBA[2]);
      }
    }
  }

  // Initialize pegs
  Pegs.clear();
  for (int idxPeg= 0; idxPeg < D.UI[PegNumber_______].I(); idxPeg++) {
    if (D.UI[PegLayout_______].I() == 0) {
      Pegs.push_back(std::array<int, 2>{Random::Val(0, nW - 1), Random::Val(0, nH - 1)});
    }
    else {
      const float angle= 2.0f * std::numbers::pi * (float)idxPeg / (float)D.UI[PegNumber_______].I();
      const int w= std::round((0.5 + 0.5 * std::cos(angle)) * (nW - 1));
      const int h= std::round((0.5 + 0.5 * std::sin(angle)) * (nH - 1));
      Pegs.push_back(std::array<int, 2>{std::min(std::max(w, 0), nW - 1), std::min(std::max(h, 0), nH - 1)});
    }
  }
  PegsCount= std::vector<int>(Pegs.size(), 0);

  // Initialize colors
  Colors.clear();
  if (D.UI[ColorsAdd_______].I() == 1) {
    Colors.push_back(Vec::Vec3<float>(1.0f, 1.0f, 1.0f));
  }
  if (D.UI[ColorsAdd_______].I() == 2) {
    Colors.push_back(Vec::Vec3<float>(1.0f, 0.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 1.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 0.0f, 1.0f));
  }
  if (D.UI[ColorsAdd_______].I() == 3) {
    Colors.push_back(Vec::Vec3<float>(0.0f, 1.0f, 1.0f));
    Colors.push_back(Vec::Vec3<float>(1.0f, 0.0f, 1.0f));
    Colors.push_back(Vec::Vec3<float>(1.0f, 1.0f, 0.0f));
  }
  if (D.UI[ColorsAdd_______].I() == 4) {
    Colors.push_back(Vec::Vec3<float>(1.0f, 1.0f, 1.0f));
    Colors.push_back(Vec::Vec3<float>(1.0f, 0.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 1.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 0.0f, 1.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 1.0f, 1.0f));
    Colors.push_back(Vec::Vec3<float>(1.0f, 0.0f, 1.0f));
    Colors.push_back(Vec::Vec3<float>(1.0f, 1.0f, 0.0f));
  }
  if (D.UI[ColorsSub_______].I() == 1) {
    Colors.push_back(Vec::Vec3<float>(-1.0f, -1.0f, -1.0f));
  }
  if (D.UI[ColorsSub_______].I() == 2) {
    Colors.push_back(Vec::Vec3<float>(-1.0f, 0.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, -1.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 0.0f, -1.0f));
  }
  if (D.UI[ColorsSub_______].I() == 3) {
    Colors.push_back(Vec::Vec3<float>(0.0f, -1.0f, -1.0f));
    Colors.push_back(Vec::Vec3<float>(-1.0f, 0.0f, -1.0f));
    Colors.push_back(Vec::Vec3<float>(-1.0f, -1.0f, 0.0f));
  }
  if (D.UI[ColorsSub_______].I() == 4) {
    Colors.push_back(Vec::Vec3<float>(-1.0f, -1.0f, -1.0f));
    Colors.push_back(Vec::Vec3<float>(-1.0f, 0.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, -1.0f, 0.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, 0.0f, -1.0f));
    Colors.push_back(Vec::Vec3<float>(0.0f, -1.0f, -1.0f));
    Colors.push_back(Vec::Vec3<float>(-1.0f, 0.0f, -1.0f));
    Colors.push_back(Vec::Vec3<float>(-1.0f, -1.0f, 0.0f));
  }

  if (D.UI[ColorsNormalize_].B())
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++)
      Colors[idxCol].normalize();

  // Initialize lines
  Lines.clear();
  Lines.resize(Colors.size(), std::vector<int>(1, 0));
}


// Handle keypress
void StringArtOptim::KeyPress(const unsigned char key) {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  (void)key;  // Disable warning unused variable
}


// Animate the project
void StringArtOptim::Animate() {
  if (!isActivProj) return;
  if (!CheckAlloc()) Allocate();
  if (!CheckRefresh()) Refresh();

  // Compute and add the new lines
  bool lineWasAdded= false;
  for (int idxStep= 0; idxStep < D.UI[StepCount_______].I(); idxStep++) {
    if (StringArtOptim::AddLineStep()) lineWasAdded= true;
    else break;
  }

  if (lineWasAdded) {
    // Compute the total match error
    float Err= 0.0f;
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        Err+= (ImCur[w][h] - ImRef[w][h]).normSquared();
      }
    }

    // Compute the total string length for each color
    std::vector<float> Len((int)Colors.size(), 0.0f);
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      for (int idxLine= 1; idxLine < (int)Lines[idxCol].size(); idxLine++) {
        Vec::Vec2<float> pos0((Pegs[Lines[idxCol][idxLine - 1]][0] + 0.5f) / (float)(nW), (Pegs[Lines[idxCol][idxLine - 1]][1] + 0.5f) / (float)(nH));
        Vec::Vec2<float> pos1((Pegs[Lines[idxCol][idxLine]][0] + 0.5f) / (float)(nW), (Pegs[Lines[idxCol][idxLine]][1] + 0.5f) / (float)(nH));
        Len[idxCol]+= (pos0 - pos1).norm();
      }
    }

    // Add to plot data
    D.Plot.resize(2);
    D.Plot[0].name= "MatchErr";
    D.Plot[1].name= "PegCounts";
    D.Plot[0].val.push_back(Err);
    D.Plot[1].val.resize(Pegs.size());
    for (int idxPeg= 0; idxPeg < (int)Pegs.size(); idxPeg++)
      D.Plot[1].val[idxPeg]= PegsCount[idxPeg];

    // Add to scatter data
    D.Scatter.resize(1);
    D.Scatter[0].name= "Lengths";
    D.Scatter[0].val.resize(Colors.size());
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      D.Scatter[0].val[idxCol]= std::array<double, 2>({(double)idxCol, (double)Len[idxCol]});
    }
  }
}


// Draw the project
void StringArtOptim::Draw() {
  if (!isActivProj) return;
  if (!isAllocated) return;
  if (!isRefreshed) return;

  // Draw the reference and generated images
  if (D.displayMode1 || D.displayMode2) {
    glPushMatrix();
    glTranslatef(0.5f, 0.0f, 0.0f);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    for (int w= 0; w < nW; w++) {
      for (int h= 0; h < nH; h++) {
        if (D.displayMode1) glColor3fv(ImCur[w][h].array());
        else if (D.displayMode2) glColor3fv(ImRef[w][h].array());
        glRectf(float(w) / float(nW), float(h) / float(nH), float(w + 1) / float(nW), float(h + 1) / float(nH));
      }
    }
    glPopMatrix();
  }

  // Draw the pegs
  if (D.displayMode3) {
    float avgPegCount= 0.0f;
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      avgPegCount+= (float)Lines[idxCol].size();
    }
    avgPegCount/= (float)Pegs.size();
    glPointSize(10.0f);
    glBegin(GL_POINTS);
    for (int idxPeg= 0; idxPeg < (int)Pegs.size(); idxPeg++) {
      float r= 0.0f, g= 0.0f, b= 0.0f;
      Colormap::RatioToJetBrightSmooth(0.5f * (float)PegsCount[idxPeg] / avgPegCount, r, g, b);
      glColor3f(r, g, b);
      glVertex3f(0.5f + 0.001f, (Pegs[idxPeg][0] + 0.5f) / float(nW), (Pegs[idxPeg][1] + 0.5f) / float(nH));
    }
    glEnd();
    glPointSize(1.0f);
  }

  // Draw the string
  if (!D.displayMode4) {
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      if (Lines[idxCol].size() >= 2) {
        glLineWidth(2.0f);
        glBegin(GL_LINE_STRIP);
        glColor3fv((0.5f * Colors[idxCol] + Vec::Vec3<float>(0.5f, 0.5f, 0.5f)).array());
        for (int idxLine= 0; idxLine < (int)Lines[idxCol].size(); idxLine++) {
          Vec::Vec3<float> pos(0.5f + 0.05f * float(idxLine) / float(Lines[idxCol].size()),
                               (Pegs[Lines[idxCol][idxLine]][0] + 0.5f) / (float)(nW),
                               (Pegs[Lines[idxCol][idxLine]][1] + 0.5f) / (float)(nH));
          glVertex3fv(pos.array());
        }
        glEnd();
        glLineWidth(1.0f);
      }
    }
  }

  // Draw the chosen string colors
  if (D.displayMode5) {
    glPushMatrix();
    glTranslatef(0.5f, 0.0f, 0.0f);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      if (Colors[idxCol].sum() > 0.0) {
        glColor3fv(Colors[idxCol].array());
        glRectf(1.0f + 0.02f, float(idxCol) * 0.1f + 0.02f, 1.0f + 0.08f, float(idxCol) * 0.1f + 0.08f);
      }
      else {
        glColor3fv((Vec::Vec3<float>(1.0f, 1.0f, 1.0f) + Colors[idxCol]).array());
        glRectf(1.1f + 0.02f, float(idxCol) * 0.1f + 0.02f, 1.1f + 0.08f, float(idxCol) * 0.1f + 0.08f);
      }
    }
    glPopMatrix();
  }
}


bool StringArtOptim::AddLineStep() {
  // Intialize the update arrays
  std::vector<int> bestPeg((int)Colors.size(), -1);
  std::vector<float> bestErr((int)Colors.size(), 0.0f);

// Find the best next peg for each color with greedy search
#pragma omp parallel for
  for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
    const int idxPeg0= Lines[idxCol].back();
    for (int idxPeg1= 0; idxPeg1 < (int)Pegs.size(); idxPeg1++) {
      std::vector<std::array<int, 2>> path= Bresenham::Line2D(Pegs[idxPeg0][0], Pegs[idxPeg0][1], Pegs[idxPeg1][0], Pegs[idxPeg1][1]);
      float chgErr= 0.0f;
      for (int idxPos= 0; idxPos < (int)path.size(); idxPos++) {
        const int w= path[idxPos][0];
        const int h= path[idxPos][1];
        Vec::Vec3<float> newVal;
        if (D.UI[BlendMode_______].B())
          newVal= (1.0f - D.UI[CoeffColor______].F()) * ImCur[w][h] + D.UI[CoeffColor______].F() * Colors[idxCol];
        else
          newVal= ImCur[w][h] + D.UI[CoeffColor______].F() * Colors[idxCol];
        for (int k= 0; k < 3; k++) newVal[k]= std::min(std::max(newVal[k], 0.0f), 1.0f);
        const Vec::Vec3<float> curErr= ImRef[w][h] - ImCur[w][h];
        const Vec::Vec3<float> newErr= ImRef[w][h] - newVal;
        chgErr+= newErr.normSquared() - curErr.normSquared();
      }
      if (bestPeg[idxCol] < 0 || bestErr[idxCol] > chgErr) {
        bestErr[idxCol]= chgErr;
        bestPeg[idxCol]= idxPeg1;
      }
    }
  }

  // Optionally keep only the best color
  if (D.UI[SingleLine______].B()) {
    int idxBestCol= 0;
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      if (bestErr[idxBestCol] > bestErr[idxCol]) {
        idxBestCol= idxCol;
      }
    }
    for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
      if (idxCol != idxBestCol) {
        bestPeg[idxCol]= -1;
      }
    }
  }

  // Update the lines
  bool lineWasAdded= false;
  for (int idxCol= 0; idxCol < (int)Colors.size(); idxCol++) {
    if (bestPeg[idxCol] >= 0 && bestErr[idxCol] < 0.0f) {
      Lines[idxCol].push_back(bestPeg[idxCol]);
      PegsCount[bestPeg[idxCol]]++;
      lineWasAdded= true;
      // Update the image
      if ((int)Lines[idxCol].size() > 1) {
        const int w0= Pegs[Lines[idxCol][Lines[idxCol].size() - 2]][0];
        const int h0= Pegs[Lines[idxCol][Lines[idxCol].size() - 2]][1];
        const int w1= Pegs[Lines[idxCol][Lines[idxCol].size() - 1]][0];
        const int h1= Pegs[Lines[idxCol][Lines[idxCol].size() - 1]][1];
        std::vector<std::array<int, 2>> path= Bresenham::Line2D(w0, h0, w1, h1);
        for (int idxPos= 0; idxPos < (int)path.size(); idxPos++) {
          const int w= path[idxPos][0];
          const int h= path[idxPos][1];
          if (D.UI[BlendMode_______].B())
            ImCur[w][h]= (1.0f - D.UI[CoeffColor______].F()) * ImCur[w][h] + D.UI[CoeffColor______].F() * Colors[idxCol];
          else
            ImCur[w][h]= ImCur[w][h] + D.UI[CoeffColor______].F() * Colors[idxCol];
          for (int k= 0; k < 3; k++) ImCur[w][h][k]= std::min(std::max(ImCur[w][h][k], 0.0f), 1.0f);
        }
      }
    }
  }

  return lineWasAdded;
}
