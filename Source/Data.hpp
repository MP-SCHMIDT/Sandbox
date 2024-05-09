#pragma once

#include <array>
#include <string>
#include <vector>


class ParamUI
{
  private:
  double val;
  bool changeFlag;

  public:
  std::string name;
  ParamUI(std::string const iName, double const iVal) {
    name= iName;
    val= iVal;
    changeFlag= true;
  }

  void Set(double const iVal) {
    changeFlag= true;
    val= iVal;
  }

  bool B() { return val > 0.0; }
  int I() { return (int)((val < 0.0) ? (val - 0.5) : (val + 0.5)); }
  float F() { return (float)val; }
  double D() { return val; }

  bool hasChanged() {
    if (changeFlag) {
      changeFlag= false;
      return true;
    }
    return false;
  }
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

  bool displayMode1= true;
  bool displayMode2= true;
  bool displayMode3= true;
  bool displayMode4= true;
  bool displayMode5= true;
  bool displayMode6= true;
  bool displayMode7= true;
  bool displayMode8= true;
  bool displayMode9= true;
  bool showAxis= true;

  std::array<double, 3> boxMin= {0.0, 0.0, 0.0};
  std::array<double, 3> boxMax= {1.0, 1.0, 1.0};

  std::array<double, 3> camDir= {0.0, 0.0, 0.0};

  std::array<double, 3> mouseNear= {0.0, 0.0, 0.0};
  std::array<double, 3> mouseFar= {0.0, 0.0, 0.0};
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
