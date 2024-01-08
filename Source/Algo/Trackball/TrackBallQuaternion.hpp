#pragma once

#include <assert.h>

#include <cmath>
#include <numbers>

#include "TrackBallVectors.hpp"

namespace tbmath {

  class quaternion
  {
public:
    enum PART
    {
      W= 0,
      X,
      Y,
      Z
    };

public:
    quaternion(float w= 0.0, float x= 0.0, float y= 0.0, float z= 0.0)
        : w(w), x(x), y(y), z(z) {}
    quaternion(const quaternion& q) { *this= q; }
    quaternion(float yaw, float pitch, float roll) { euler(yaw, pitch, roll); }
    quaternion(const vec3& v) { set(0.0, v[0], v[1], v[2]); }
    quaternion& set(float a, float b, float c, float d) {
      w= a;
      x= b;
      y= c;
      z= d;
      return *this;
    }
    quaternion(float angle, const vec3& axis) {
      float d= axis.length();
      assert(d != 0.0 && "quaternion Axis is Zero");
      angle*= (std::numbers::pi / 180.0);
      float s= sin(angle * 0.5) / d;
      set(cos(angle * 0.5), axis[0] * s, axis[1] * s, axis[2] * s);
    }
    quaternion& euler(float yaw, float pitch, float roll) {
      yaw*= (std::numbers::pi / 180.0) / 2.0;
      pitch*= (std::numbers::pi / 180.0) / 2.0;
      roll*= (std::numbers::pi / 180.0) / 2.0;
      float c_yaw= cos(yaw);
      float s_yaw= sin(yaw);
      float c_pitch= cos(pitch);
      float s_pitch= sin(pitch);
      float c_roll= cos(roll);
      float s_roll= sin(roll);
      set(
          c_roll * c_pitch * c_yaw + s_roll * s_pitch * s_yaw,
          c_roll * s_pitch * c_yaw + s_roll * c_pitch * s_yaw,
          c_roll * c_pitch * s_yaw - s_roll * s_pitch * c_yaw,
          s_roll * c_pitch * c_yaw - c_roll * s_pitch * s_yaw);
      return *this;
    }
    inline float dot(const quaternion& q) { return (w * q.w + x * q.x + y * q.y + z * q.z); }
    inline quaternion& normalize(float len= 1.0) {
      float l= length();
      return (l != 0.0) ? (*this*= (len / l)) : (*this);
    }
    inline float length() const { return sqrt(length2()); }
    inline float length2() const { return w * w + x * x + y * y + z * z; }
    inline vec3 getvec3() const { return vec3(x, y, z); }
    inline vec3& rot(vec3& v) {
      v= ((*this * quaternion(v)) * quaternion(w, -x, -y, -z))
             .getvec3();
      v*= 1.0 / length2();
      return v;
    }
    inline float operator[](int idx) const {
      switch (idx) {
        case W:
          return w;
        case X:
          return x;
        case Y:
          return y;
        case Z:
          return z;
      }
      return 0.0;
    }
    inline quaternion& operator=(const quaternion& q) {
      w= q.w;
      x= q.x;
      y= q.y;
      z= q.z;
      return *this;
    }
    inline bool operator==(const quaternion& q) { return (w == q.w) && (x == q.x) && (y == q.y) && (z == q.z); }
    inline quaternion& operator*=(const quaternion& q) { return *this= quaternion(*this)*= q; }
    inline quaternion& operator*=(float s) {
      w*= s;
      x*= s;
      y*= s;
      z*= s;
      return *this;
    }
    inline quaternion& operator+=(const quaternion& q) {
      w+= q.w;
      x+= q.x;
      y+= q.y;
      z+= q.z;
      return *this;
    }
    inline quaternion& operator-=(const quaternion& q) {
      w-= q.w;
      x-= q.x;
      y-= q.y;
      z-= q.z;
      return *this;
    }
    friend quaternion operator+(const quaternion& a, const quaternion& b) { return quaternion(a)+= b; }
    friend quaternion operator-(const quaternion& a, const quaternion& b) { return quaternion(a)-= b; }
    friend quaternion operator*(const quaternion& a, const quaternion& b) {
      return quaternion(
          a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
          a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
          a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
          a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x);
    }
    friend quaternion operator*(const quaternion& q, float s) { return quaternion(q)*= s; }
    friend quaternion operator*(float s, const quaternion& q) { return quaternion(q)*= s; }
    friend quaternion operator*(const vec3& v, const quaternion& q) { return quaternion(0, v)*= q; }
    friend quaternion operator*(const quaternion& q, const vec3& v) { return q * quaternion(0, v); }

private:
    float w, x, y, z;
  };

}  // namespace tbmath
