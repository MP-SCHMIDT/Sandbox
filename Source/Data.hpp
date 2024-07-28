#pragma once

#include <array>
#include <string>
#include <vector>


class ParamUI
{
  public:
  double val;
  bool changeFlag;

  std::string name;
  ParamUI(std::string const iName, double const iVal) {
    name= iName;
    val= iVal;
    changeFlag= true;
  }

  inline void Set(double const iVal) {
    changeFlag= true;
    val= iVal;
  }

  inline bool hasChanged(const bool iResetFlag= true) {
    if (changeFlag) {
      if (iResetFlag) changeFlag= false;
      return true;
    }
    return false;
  }

  inline bool B() { return val >= 0.5; }
  inline int I() { return (int)((val < 0.0) ? (val - 0.5) : (val + 0.5)); }
  inline float F() { return (float)val; }
  inline double D() { return val; }
};


class PlotUI
{
  public:
  std::vector<double> val;
  std::string name;
  bool isLog;
  bool isSameRange;
  bool isSymmetric;
  bool showPoints;

  PlotUI() {
    val.clear();
    name= "<name>";
    isLog= false;
    isSameRange= false;
    isSymmetric= false;
    showPoints= false;
  }
};


class ScatterUI
{
  public:
  std::vector<std::array<double, 2>> val;
  std::string name;

  ScatterUI() {
    val.clear();
    name= "<name>";
  }
};


class Data
{
  public:
  bool playAnimation= false;
  bool stepAnimation= false;

  std::array<bool, 10> displayMode= {true, true, true, true, true, true, true, true, true, true};
  std::array<std::string, 10> displayModeLabel;

  bool showAxis= true;
  std::array<double, 3> boxMin= {0.0, 0.0, 0.0};
  std::array<double, 3> boxMax= {1.0, 1.0, 1.0};

  std::array<double, 3> camDir= {0.0, 0.0, 0.0};

  bool keyIsShift= false;
  bool keyIsCtrl= false;
  bool keyIsAlt= false;
  unsigned char keyLetterUpperCase= 0;
  unsigned char mouseMiddleButtonState= 0;
  std::array<double, 3> mouseBeg= {0.0, 0.0, 0.0};
  std::array<double, 3> mouseEnd= {0.0, 0.0, 0.0};
  std::array<double, 3> mouseProjX= {0.0, 0.0, 0.0};
  std::array<double, 3> mouseProjY= {0.0, 0.0, 0.0};
  std::array<double, 3> mouseProjZ= {0.0, 0.0, 0.0};

  int idxFirstParamPageUI= 0;  // Index of currently displayed parameter page
  int idxParamUI= 0;           // Index of currently selected parameter
  int idxCursorUI= 0;          // Index of currently selected character in the parameter value

  std::vector<ParamUI> UI;
  std::vector<PlotUI> Plot;
  std::vector<ScatterUI> Scatter;
  std::vector<std::string> Status;
};
