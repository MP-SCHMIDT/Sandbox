#pragma once


class DrawShape
{
  public:
  static void DrawBoxPosPos(const float begX, const float begY, const float begZ,
                            const float endX, const float endY, const float endZ, bool const isSolid);

  static void DrawBoxPosSiz(const float begX, const float begY, const float begZ,
                            const float sizX, const float sizY, const float sizZ, bool const isSolid);
};
