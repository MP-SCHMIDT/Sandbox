#pragma once

// Standard lib
#include <cmath>
#include <concepts>


namespace Vec {

  template <std::floating_point element_type>
  class Vec2
  {
public:
    // Constructors
    Vec2(element_type a= 0.0, element_type b= 0.0) {
      this->x[0]= a;
      this->x[1]= b;
    }
    Vec2(const Vec2& v) { *this= v; }

    // Methods
    inline void set(element_type a, element_type b) {
      x[0]= a;
      x[1]= b;
    }
    inline element_type dot(const Vec2& v) const { return x[0] * v.x[0] + x[1] * v.x[1]; }
    inline Vec2 coeffMul(const Vec2& v) const { return Vec2(x[0] * v.x[0], x[1] * v.x[1]); }
    inline Vec2 coeffDiv(const Vec2& v) const { return Vec2(x[0] / v.x[0], x[1] / v.x[1]); }
    inline element_type normSquared() const { return x[0] * x[0] + x[1] * x[1]; }
    inline element_type norm() const { return sqrt(normSquared()); }
    inline Vec2& normalize(element_type len= 1.0) {
      element_type l= norm();
      return (l != 0.0) ? (*this*= (len / l)) : (*this);
    }
    inline Vec2 normalized(element_type len= 1.0) {
      element_type l= norm();
      return (l != 0.0) ? (*this * len / l) : (*this);
    }
    inline Vec2 abs() {
      return Vec2(std::abs(x[0]), std::abs(x[1]));
    }
    inline element_type sum() {
      return x[0] + x[1];
    }
    inline element_type max() const { return std::max(x[0], x[1]); }
    inline element_type min() const { return std::min(x[0], x[1]); }

    // Operators
    inline element_type operator[](int idx) const { return x[idx]; }
    inline element_type& operator[](int idx) { return x[idx]; }
    inline operator const element_type*(void) const { return x; }
    inline Vec2& operator=(const Vec2& v) {
      x[0]= v.x[0];
      x[1]= v.x[1];
      return *this;
    }
    inline bool operator==(const Vec2& v) { return ((x[0] == v.x[0]) && (x[1] == v.x[1])); }
    inline Vec2& operator+=(const Vec2& v) {
      x[0]+= v.x[0];
      x[1]+= v.x[1];
      return *this;
    }
    inline Vec2& operator-=(const Vec2& v) {
      x[0]-= v.x[0];
      x[1]-= v.x[1];
      return *this;
    }
    inline Vec2& operator*=(element_type f) {
      x[0]*= f;
      x[1]*= f;
      return *this;
    }
    inline Vec2& operator/=(element_type f) {
      x[0]/= f;
      x[1]/= f;
      return *this;
    }
    friend Vec2 operator+(const Vec2& w, const Vec2& v) { return Vec2(w)+= v; }
    friend Vec2 operator-(const Vec2& w, const Vec2& v) { return Vec2(w)-= v; }
    friend Vec2 operator*(const Vec2& v, element_type f) { return Vec2(v)*= f; }
    friend Vec2 operator*(element_type f, const Vec2& v) { return Vec2(v)*= f; }
    friend Vec2 operator/(const Vec2& v, element_type f) { return Vec2(v)/= f; }

    inline element_type* array() { return x; }

private:
    element_type x[2];
  };


  template <std::floating_point element_type>
  class Vec3
  {
public:
    // Constructors
    Vec3(element_type a= 0.0, element_type b= 0.0, element_type c= 0.0) {
      this->x[0]= a;
      this->x[1]= b;
      this->x[2]= c;
    }
    Vec3(const Vec3& v) { *this= v; }

    // Methods
    inline void set(element_type a, element_type b, element_type c) {
      x[0]= a;
      x[1]= b;
      x[2]= c;
    }
    inline element_type dot(const Vec3& v) const { return x[0] * v.x[0] + x[1] * v.x[1] + x[2] * v.x[2]; }
    inline Vec3 cross(const Vec3& v) const {
      return Vec3(
          x[1] * v.x[2] - x[2] * v.x[1],
          x[2] * v.x[0] - x[0] * v.x[2],
          x[0] * v.x[1] - x[1] * v.x[0]);
    }
    inline Vec3 coeffMul(const Vec3& v) const { return Vec3(x[0] * v.x[0], x[1] * v.x[1], x[2] * v.x[2]); }
    inline Vec3 coeffDiv(const Vec3& v) const { return Vec3(x[0] / v.x[0], x[1] / v.x[1], x[2] / v.x[2]); }
    inline element_type normSquared() const { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; }
    inline element_type norm() const { return sqrt(normSquared()); }
    inline Vec3& normalize(element_type len= 1.0) {
      element_type l= norm();
      return (l != 0.0) ? (*this*= (len / l)) : (*this);
    }
    inline Vec3 normalized(element_type len= 1.0) {
      element_type l= norm();
      return (l != 0.0) ? (*this * len / l) : (*this);
    }
    inline Vec3 abs() {
      return Vec3(std::abs(x[0]), std::abs(x[1]), std::abs(x[2]));
    }
    inline element_type sum() {
      return x[0] + x[1] + x[2];
    }
    inline element_type max() const { return std::max(std::max(x[0], x[1]), x[2]); }
    inline element_type min() const { return std::min(std::min(x[0], x[1]), x[2]); }

    // Operators
    inline element_type operator[](int idx) const { return x[idx]; }
    inline element_type& operator[](int idx) { return x[idx]; }
    inline operator const element_type*(void) const { return x; }
    inline Vec3& operator=(const Vec3& v) {
      x[0]= v.x[0];
      x[1]= v.x[1];
      x[2]= v.x[2];
      return *this;
    }
    inline bool operator==(const Vec3& v) { return ((x[0] == v.x[0]) && (x[1] == v.x[1]) && (x[2] == v.x[2])); }
    inline Vec3& operator+=(const Vec3& v) {
      x[0]+= v.x[0];
      x[1]+= v.x[1];
      x[2]+= v.x[2];
      return *this;
    }
    inline Vec3& operator-=(const Vec3& v) {
      x[0]-= v.x[0];
      x[1]-= v.x[1];
      x[2]-= v.x[2];
      return *this;
    }
    inline Vec3& operator*=(element_type f) {
      x[0]*= f;
      x[1]*= f;
      x[2]*= f;
      return *this;
    }
    inline Vec3& operator/=(element_type f) {
      x[0]/= f;
      x[1]/= f;
      x[2]/= f;
      return *this;
    }
    friend Vec3 operator+(const Vec3& w, const Vec3& v) { return Vec3(w)+= v; }
    friend Vec3 operator-(const Vec3& w, const Vec3& v) { return Vec3(w)-= v; }
    friend Vec3 operator*(const Vec3& v, element_type f) { return Vec3(v)*= f; }
    friend Vec3 operator*(element_type f, const Vec3& v) { return Vec3(v)*= f; }
    friend Vec3 operator/(const Vec3& v, element_type f) { return Vec3(v)/= f; }

    inline element_type* array() { return x; }

private:
    element_type x[3];
  };


  template <std::floating_point element_type>
  class Vec4
  {
public:
    // Constructors
    Vec4(element_type a= 0.0, element_type b= 0.0, element_type c= 0.0, element_type d= 0.0) {
      this->x[0]= a;
      this->x[1]= b;
      this->x[2]= c;
      this->x[3]= d;
    }
    Vec4(const Vec4& v) { *this= v; }

    // Methods
    inline void set(element_type a, element_type b, element_type c, element_type d) {
      x[0]= a;
      x[1]= b;
      x[2]= c;
      x[3]= d;
    }
    inline element_type dot(const Vec4& v) const { return x[0] * v.x[0] + x[1] * v.x[1] + x[2] * v.x[2] + x[3] * v.x[3]; }
    inline Vec4 coeffMul(const Vec4& v) const { return Vec4(x[0] * v.x[0], x[1] * v.x[1], x[2] * v.x[2], x[3] * v.x[3]); }
    inline Vec4 coeffDiv(const Vec4& v) const { return Vec4(x[0] / v.x[0], x[1] / v.x[1], x[2] / v.x[2], x[3] / v.x[3]); }
    inline element_type normSquared() const { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; }
    inline element_type norm() const { return sqrt(normSquared()); }
    inline Vec4& normalize(element_type len= 1.0) {
      element_type l= norm();
      return (l != 0.0) ? (*this*= (len / l)) : (*this);
    }
    inline Vec4 normalized(element_type len= 1.0) {
      element_type l= norm();
      return (l != 0.0) ? (*this * len / l) : (*this);
    }
    inline Vec4 abs() {
      return Vec4(std::abs(x[0]), std::abs(x[1]), std::abs(x[2]), std::abs(x[3]));
    }
    inline element_type sum() {
      return x[0] + x[1] + x[2] + x[3];
    }
    inline element_type max() const { return std::max(std::max(std::max(x[0], x[1]), x[2]), x[3]); }
    inline element_type min() const { return std::min(std::min(std::min(x[0], x[1]), x[2]), x[3]); }

    // Operators
    inline element_type operator[](int idx) const { return x[idx]; }
    inline element_type& operator[](int idx) { return x[idx]; }
    inline operator const element_type*(void) const { return x; }
    inline Vec4& operator=(const Vec4& v) {
      x[0]= v.x[0];
      x[1]= v.x[1];
      x[2]= v.x[2];
      x[3]= v.x[3];
      return *this;
    }
    inline bool operator==(const Vec4& v) { return ((x[0] == v.x[0]) && (x[1] == v.x[1]) && (x[2] == v.x[2]) && (x[3] == v.x[3])); }
    inline Vec4& operator+=(const Vec4& v) {
      x[0]+= v.x[0];
      x[1]+= v.x[1];
      x[2]+= v.x[2];
      x[3]+= v.x[3];
      return *this;
    }
    inline Vec4& operator-=(const Vec4& v) {
      x[0]-= v.x[0];
      x[1]-= v.x[1];
      x[2]-= v.x[2];
      x[3]-= v.x[3];
      return *this;
    }
    inline Vec4& operator*=(element_type f) {
      x[0]*= f;
      x[1]*= f;
      x[2]*= f;
      x[3]*= f;
      return *this;
    }
    inline Vec4& operator/=(element_type f) {
      x[0]/= f;
      x[1]/= f;
      x[2]/= f;
      x[3]/= f;
      return *this;
    }
    friend Vec4 operator+(const Vec4& w, const Vec4& v) { return Vec4(w)+= v; }
    friend Vec4 operator-(const Vec4& w, const Vec4& v) { return Vec4(w)-= v; }
    friend Vec4 operator*(const Vec4& v, element_type f) { return Vec4(v)*= f; }
    friend Vec4 operator*(element_type f, const Vec4& v) { return Vec4(v)*= f; }
    friend Vec4 operator/(const Vec4& v, element_type f) { return Vec4(v)/= f; }

    inline element_type* array() { return x; }

private:
    element_type x[4];
  };

}  // namespace Vec
