#pragma once

#include "TrackBallMatrix.hpp"
#include "TrackBallQuaternion.hpp"
#include "TrackBallVectors.hpp"

namespace tb {

  // TrackBall implementation, inspired by Ken Shoemake
  class TrackBall
  {
private:
    float width;               // width of window
    float height;              // height of window
    tbmath::vec3 center;       // center of the virtual trackball
    double radius;             // radius of the virtual trackball
    bool drag;                 // flag: true=dragging, false=nothing
    tbmath::vec3 v_cur;        // current mouse position
    tbmath::vec3 v_down;       // mouse position at the beginning of dragging
    tbmath::quaternion q_cur;  // current rotation
    tbmath::quaternion q_end;  // rotation after the dragging
    tbmath::mat4 mat_cur;      // current rotation matrix
private:
    /// maps the specified mouse position to the sphere defined
    /// with center and radius. the resulting vector lies on the
    /// surface of the sphere.
    inline tbmath::vec3 map_sphere(const tbmath::vec3& mouse,
                                   const tbmath::vec3& center, double r) const {
      tbmath::vec3 bm;

      bm= (1.0 / r) * (mouse - center);
      double mag= bm.length2();
      if (mag > 1.0) {
        bm.normalize();
        bm[2]= 0.0;
      }
      else {
        bm[2]= std::sqrt(1.0 - mag);
      }
      return bm;
    }

private:
    /// calculates and returns the quaternion which describes the
    /// arc between the two specified positions on the surface
    /// of a sphere (virtual trackball).
    inline tbmath::quaternion from_ball_points(
        const tbmath::vec3& from, const tbmath::vec3& to) const {
      return tbmath::quaternion(
          from[0] * to[0] + from[1] * to[1] + from[2] * to[2], from[1] * to[2] - from[2] * to[1], from[2] * to[0] - from[0] * to[2], from[0] * to[1] - from[1] * to[0]);
    }

private:
    /// maps the specified quaternion to the 4x4 matrix
    /// which represents a rotation. works for right-handed
    /// coordinate systems.
    inline void to_matrix(tbmath::mat4& m,
                          const tbmath::quaternion& q) const {
      double l= q.length2();
      double s= (l > 0.0) ? (2.0 / l) : 0.0;

      double xs= q[1] * s;
      double ys= q[2] * s;
      double zs= q[3] * s;

      double wx= q[0] * xs;
      double wy= q[0] * ys;
      double wz= q[0] * zs;

      double xx= q[1] * xs;
      double xy= q[1] * ys;
      double xz= q[1] * zs;

      double yy= q[2] * ys;
      double yz= q[2] * zs;

      double zz= q[3] * zs;

      m[0]= 1.0 - (yy + zz);
      m[1]= xy + wz;
      m[2]= xz - wy;
      m[3]= 0.0;

      m[4]= xy - wz;
      m[5]= 1.0 - (xx + zz);
      m[6]= yz + wx;
      m[7]= 0.0;

      m[8]= xz + wy;
      m[9]= yz - wx;
      m[10]= 1.0 - (xx + yy);
      m[11]= 0.0;

      m[12]= 0.0;
      m[13]= 0.0;
      m[14]= 0.0;
      m[15]= 1.0;
    }

private:
    /// update of the rotation matrix 'mat_cur', using the position
    /// of the beginning of the dragging and the current position.
    /// both coordinates mapped to the surface of the virtual trackball.
    inline void update(void) {
      tbmath::vec3 v_from= map_sphere(v_down, center, radius);
      tbmath::vec3 v_to= map_sphere(v_cur, center, radius);
      if (drag) q_cur= from_ball_points(v_from, v_to) * q_end;
      to_matrix(mat_cur, q_cur);
    }

public:
    TrackBall(void)
        : center(0.0, 0.0, 0.0), radius(1.0), drag(false), v_cur(0.0, 0.0, 0.0), v_down(0.0, 0.0, 0.0), q_cur(1.0, 0.0, 0.0, 0.0), q_end(1.0, 0.0, 0.0, 0.0) {}

public:
    /// sets the window size.
    inline void set_win_size(float width, float height) {
      this->width= width;
      this->height= height;
    }

public:
    /// moves the center of the virtual trackball to the specified
    /// position.
    inline void place(const tbmath::vec3& center) { this->center= center; }

public:
    /// sets the radius of the virtual trackball.
    inline void set_r(double r) { this->radius= r; }

public:
    /// sets the current position and calculates the current
    /// rotation matrix.
    inline void set_cur(int x, int y) {
      v_cur[0]= 2.0 * ((double)x / width) - 1.0;
      v_cur[1]= 2.0 * ((double)(height - y) / height) - 1.0;
      v_cur[2]= 0.0;
      update();  // calculate rotation matrix
    }

public:
    /// indicates the beginning of the dragging.
    inline void begin_drag(void) {
      drag= true;     // start dragging
      v_down= v_cur;  // remember start position
    }

public:
    /// marks the end of the dragging.
    inline void end_drag(void) {
      drag= false;   // stop dragging
      q_end= q_cur;  // remember rotation
    }

public:
    /// returns the rotation matrix to be used directly
    /// for the OpenGL command 'glMultMatrixf'.
    inline const float* get(void) const { return mat_cur.get(); }

    inline const tbmath::mat4 get_mat(void) const {
      return mat_cur;
    }
  };

}  // namespace tb
