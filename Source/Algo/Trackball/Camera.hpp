#pragma once

#include "TrackBall.hpp"

//  Represents a camera. Enables the user to move the camera and at the end
//  get the according view transformation matrix.
class Camera {
 private:
  // The arcball used for rotation with mouse movement
  tb::TrackBall arcBall_;

  // The current transformation matrix
  tbmath::mat4 viewMatrix_;

  // Vector describing the direction from the eye to the center.
  // I store this vector rather than the eye position because I want
  // to be able to move the camera with the center (needed for panning)
  tbmath::vec3 eyeDir_;

  // The point in space where the camera looks at
  tbmath::vec3 center_;

  // A factor multiplied with zoom deltas to compensate for different
  // scroll wheel speeds.
  float zoomSpeed_;

  // The current width and height of the window. Needed for accurate behaivour
  // according to the mouse movements
  float windowWidth_;
  float windowHeight_;

  // The last and current mouse position (x,y) in window space. Needed for
  // getting the relative mouse movement for panning and rotation
  tbmath::vec2 lastMousePos_;
  tbmath::vec2 currentMousePos_;

  // If true the user currently is in panning or zooming mode
  enum Mode { IDLE,
              PANNING,
              ZOOMING };
  Mode currentMode_;

 public:
  // ========================= CAMERA SETUP ================================

  // Creates a new camera with eye position (0,0,-1) looking at (0,0,0)
  // the up vector is always (0,1,0)
  Camera();

  // Creates a new camera with specified eye and center position
  // the up vector is always (0,1,0)
  Camera(float eyex, float eyey, float eyez, float centerx, float centery, float centerz);

  void setCenter(float x, float y, float z);
  void setEye(float x, float y, float z);

  void setCurrentMousePos(float x, float y);
  void setWindowSize(float w, float h);

  // ========================= CAMERA MOVEMENT =============================

  // Begins and ends the rotation state. Call this method when user clicks
  void beginRotate();
  void endRotate();

  // Begins and ends the panning state. Call this method when user clicks
  // for example the right button.
  void beginPan();
  void endPan();

  // Begins and ends the zooming state. Call this method when user clicks
  // for example the middle button.
  void beginZoom();
  void endZoom();

  // Directly moves the center by deltax in the x direction and deltay in
  // the y direction of the current viewing plane. Use this method for panning
  // with the keyboard. For panning with mouse interaction rather use the
  // combination: setCurrentMousePos() -> beginPan() -> setCurrentMousePos()
  //  -> endPan()
  void pan(float deltax, float deltay);

  // Moves the camera away from or towards the center. If you use some mouse
  // movement (not the scroll wheel) to scroll, rather use beginZoom() and
  // endZoom()
  void zoom(float deltaz);

  // ========================= MATRIX GETTER ===============================

  // Get the current camera matrix in the right format for opengl. Load this
  // matrix into the OpenGL context before drawing to get the right camera
  // view.
  const float* getViewMatrix();

 private:
  void updatePanningPosition();

  // Updates the current view matrix
  void updateViewMatrix();

  // Gets the up and right vector from the current view matrix
  tbmath::vec3 currentUp();
  tbmath::vec3 currentRight();

  void initConstants();
};
