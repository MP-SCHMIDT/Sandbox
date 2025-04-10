#pragma once


class Colormap
{
  public:
  static float ApplyWave(const float iVal, const float iMultip);

  static void RatioToGrayscale(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToPlasma(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToViridis(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToBlackBody(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToBlueToRed(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToGreenToRed(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToRainbow(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToRedGreenBlueRed(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToJet(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToJetSmooth(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToJetBright(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToJetBrightSmooth(const float &iRatio, float &oR, float &oG, float &oB);
  static void RatioToTurbo(const float &iRatio, float &oR, float &oG, float &oB);
};
