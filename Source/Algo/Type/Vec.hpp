#pragma once

// Standard lib
#include <array>
#include <cmath>
#include <concepts>
#include <algorithm>


namespace Vec {

  template <std::floating_point element_type>
  class Vec2
  {
    public:
    // Constructors
    Vec2() { this->x[0]= 0.0; this->x[1]= 0.0; }
    Vec2(element_type a) { this->x[0]= a; this->x[1]= a; }
    Vec2(element_type a, element_type b) { this->x[0]= a; this->x[1]= b; }
    Vec2(const std::array<element_type, 2>& v) { this->x[0]= v[0]; this->x[1]= v[1]; }
    Vec2(const Vec2& v) { *this= v; }

    // Methods
    inline void set(element_type a, element_type b) { x[0]= a; x[1]= b; }
    inline element_type dot(const Vec2& v) const { return x[0] * v.x[0] + x[1] * v.x[1]; }
    inline Vec2 cwiseMul(const Vec2& v) const { return Vec2(x[0] * v.x[0], x[1] * v.x[1]); }
    inline Vec2 cwiseDiv(const Vec2& v) const { return Vec2(x[0] / v.x[0], x[1] / v.x[1]); }
    inline Vec2 cwiseMax(const Vec2& v) const { return Vec2(std::max(x[0], v.x[0]), std::max(x[1], v.x[1])); }
    inline Vec2 cwiseMin(const Vec2& v) const { return Vec2(std::min(x[0], v.x[0]), std::min(x[1], v.x[1])); }
    inline element_type normSquared() const { return x[0] * x[0] + x[1] * x[1]; }
    inline element_type norm() const { return std::sqrt(normSquared()); }
    inline void normalize(const element_type len= 1.0) { *this*= (len / norm()); }
    inline Vec2 normalized(const element_type len= 1.0) const { return Vec2(*this * len / norm()); }
    inline Vec2 abs() const { return Vec2(std::abs(x[0]), std::abs(x[1])); }
    inline element_type sum() const { return x[0] + x[1]; }
    inline element_type maxCoeff() const { return std::max(x[0], x[1]); }
    inline element_type minCoeff() const { return std::min(x[0], x[1]); }

    // Operators
    inline element_type operator[](int idx) const { return x[idx]; }
    inline element_type& operator[](int idx) { return x[idx]; }
    inline operator const element_type*(void) const { return x; }
    inline Vec2& operator=(const Vec2& v) { x[0]= v.x[0]; x[1]= v.x[1]; return *this; }
    inline bool operator==(const Vec2& v) const { return ((x[0] == v.x[0]) && (x[1] == v.x[1])); }
    inline Vec2& operator+=(const Vec2& v) { x[0]+= v.x[0]; x[1]+= v.x[1]; return *this; }
    inline Vec2& operator-=(const Vec2& v) { x[0]-= v.x[0]; x[1]-= v.x[1]; return *this; }
    inline Vec2& operator*=(const element_type f) { x[0]*= f; x[1]*= f; return *this; }
    inline Vec2& operator/=(const element_type f) { x[0]/= f; x[1]/= f; return *this; }
    friend Vec2 operator+(Vec2 v, const Vec2& w) noexcept { return v+= w; }
    friend Vec2 operator-(Vec2 v, const Vec2& w) noexcept { return v-= w; }
    friend Vec2 operator*(Vec2 v, const element_type f) noexcept { return v*= f; }
    friend Vec2 operator*(const element_type f, Vec2 v) noexcept { return v*= f; }
    friend Vec2 operator/(Vec2 v, const element_type f) noexcept { return v/= f; }

    inline element_type* array() { return x; }

    private:
    element_type x[2];
  };


  template <std::floating_point element_type>
  class Vec3
  {
    public:
    // Constructors
    Vec3() { this->x[0]= 0.0; this->x[1]= 0.0; this->x[2]= 0.0; }
    Vec3(element_type a) { this->x[0]= a; this->x[1]= a; this->x[2]= a; }
    Vec3(element_type a, element_type b, element_type c) { this->x[0]= a; this->x[1]= b; this->x[2]= c; }
    Vec3(const std::array<element_type, 3>& v) { this->x[0]= v[0]; this->x[1]= v[1]; this->x[2]= v[2]; }
    Vec3(const Vec3& v) { *this= v; }
  
    // Methods
    inline void set(element_type a, element_type b, element_type c) { x[0]= a; x[1]= b; x[2]= c; }
    inline element_type dot(const Vec3& v) const { return x[0] * v.x[0] + x[1] * v.x[1] + x[2] * v.x[2]; }
    inline Vec3 cwiseMul(const Vec3& v) const { return Vec3(x[0] * v.x[0], x[1] * v.x[1], x[2] * v.x[2]); }
    inline Vec3 cwiseDiv(const Vec3& v) const { return Vec3(x[0] / v.x[0], x[1] / v.x[1], x[2] / v.x[2]); }
    inline Vec3 cwiseMax(const Vec3& v) const { return Vec3(std::max(x[0], v.x[0]), std::max(x[1], v.x[1]), std::max(x[2], v.x[2])); }
    inline Vec3 cwiseMin(const Vec3& v) const { return Vec3(std::min(x[0], v.x[0]), std::min(x[1], v.x[1]), std::min(x[2], v.x[2])); }
    inline element_type normSquared() const { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; }
    inline element_type norm() const { return std::sqrt(normSquared()); }
    inline void normalize(const element_type len= 1.0) { *this*= (len / norm()); }
    inline Vec3 normalized(const element_type len= 1.0) const { return Vec3(*this * len / norm()); }
    inline Vec3 abs() const { return Vec3(std::abs(x[0]), std::abs(x[1]), std::abs(x[2])); }
    inline element_type sum() const { return x[0] + x[1] + x[2]; }
    inline element_type maxCoeff() const { return std::max({x[0], x[1], x[2]}); }
    inline element_type minCoeff() const { return std::min({x[0], x[1], x[2]}); }

    // Methods specific to 3D vectors
    inline Vec3 cross(const Vec3& v) const {
      return Vec3(
          x[1] * v.x[2] - x[2] * v.x[1],
          x[2] * v.x[0] - x[0] * v.x[2],
          x[0] * v.x[1] - x[1] * v.x[0]);
    }
    inline void computeBasis(Vec3& oDir2, Vec3& oDir3) {
      oDir2= this->cross(Vec3<element_type>(1.0, 0.0, 0.0));
      if (oDir2.normSquared() == 0.0)
        oDir2= this->cross(Vec3<element_type>(0.0, 1.0, 0.0));
      oDir2.normalize();
      oDir3= this->cross(oDir2);
      oDir3.normalize();
    }
    inline void computeBasisStable(Vec3& oDir2, Vec3& oDir3) {
      oDir2= this->cross(Vec3<element_type>(1.0, 0.0, 0.0));
      Vec3<element_type> const tmpVec= this->cross(Vec3<element_type>(0.0, 1.0, 0.0));
      if (oDir2.normSquared() < tmpVec.normSquared())
        oDir2= tmpVec;
      oDir2.normalize();
      oDir3= this->cross(oDir2);
      oDir3.normalize();
    }
  
    // Operators
    inline element_type operator[](int idx) const { return x[idx]; }
    inline element_type& operator[](int idx) { return x[idx]; }
    inline operator const element_type*(void) const { return x; }
    inline Vec3& operator=(const Vec3& v) { x[0]= v.x[0]; x[1]= v.x[1]; x[2]= v.x[2]; return *this; }
    inline bool operator==(const Vec3& v) const { return ((x[0] == v.x[0]) && (x[1] == v.x[1]) && (x[2] == v.x[2])); }
    inline Vec3& operator+=(const Vec3& v) { x[0]+= v.x[0]; x[1]+= v.x[1]; x[2]+= v.x[2]; return *this; }
    inline Vec3& operator-=(const Vec3& v) { x[0]-= v.x[0]; x[1]-= v.x[1]; x[2]-= v.x[2]; return *this; }
    inline Vec3& operator*=(const element_type f) { x[0]*= f; x[1]*= f; x[2]*= f; return *this; }
    inline Vec3& operator/=(const element_type f) { x[0]/= f; x[1]/= f; x[2]/= f; return *this; }
    friend Vec3 operator+(Vec3 v, const Vec3& w) noexcept { return v+= w; }
    friend Vec3 operator-(Vec3 v, const Vec3& w) noexcept { return v-= w; }
    friend Vec3 operator*(Vec3 v, const element_type f) noexcept { return v*= f; }
    friend Vec3 operator*(const element_type f, Vec3 v) noexcept { return v*= f; }
    friend Vec3 operator/(Vec3 v, const element_type f) noexcept { return v/= f; }

    inline element_type* array() { return x; }

    private:
    element_type x[3];
  };


  template <std::floating_point element_type>
  class Vec4
  {
    public:
    // Constructors
    Vec4() { this->x[0]= 0.0; this->x[1]= 0.0; this->x[2]= 0.0; this->x[3]= 0.0; }
    Vec4(element_type a) { this->x[0]= a; this->x[1]= a; this->x[2]= a; this->x[3]= a; }
    Vec4(element_type a, element_type b, element_type c, element_type d) { this->x[0]= a; this->x[1]= b; this->x[2]= c; this->x[3]= d; }
    Vec4(const std::array<element_type, 4>& v) { this->x[0]= v[0]; this->x[1]= v[1]; this->x[2]= v[2]; this->x[3]= v[3]; }
    Vec4(const Vec4& v) { *this= v; }

    // Methods
    inline void set(element_type a, element_type b, element_type c, element_type d) { x[0]= a; x[1]= b; x[2]= c; x[3]= d; }
    inline element_type dot(const Vec4& v) const { return x[0] * v.x[0] + x[1] * v.x[1] + x[2] * v.x[2] + x[3] * v.x[3]; }
    inline Vec4 cwiseMul(const Vec4& v) const { return Vec4(x[0] * v.x[0], x[1] * v.x[1], x[2] * v.x[2], x[3] * v.x[3]); }
    inline Vec4 cwiseDiv(const Vec4& v) const { return Vec4(x[0] / v.x[0], x[1] / v.x[1], x[2] / v.x[2], x[3] / v.x[3]); }
    inline Vec4 cwiseMax(const Vec4& v) const { return Vec4(std::max(x[0], v.x[0]), std::max(x[1], v.x[1]), std::max(x[2], v.x[2]), std::max(x[3], v.x[3])); }
    inline Vec4 cwiseMin(const Vec4& v) const { return Vec4(std::min(x[0], v.x[0]), std::min(x[1], v.x[1]), std::min(x[2], v.x[2]), std::min(x[3], v.x[3])); }
    inline element_type normSquared() const { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; }
    inline element_type norm() const { return std::sqrt(normSquared()); }
    inline void normalize(const element_type len= 1.0) { *this*= (len / norm()); }
    inline Vec4 normalized(const element_type len= 1.0) const { return Vec4(*this * len / norm()); }
    inline Vec4 abs() const { return Vec4(std::abs(x[0]), std::abs(x[1]), std::abs(x[2]), std::abs(x[3])); }
    inline element_type sum() const { return x[0] + x[1] + x[2] + x[3]; }
    inline element_type maxCoeff() const { return std::max({x[0], x[1], x[2], x[3]}); }
    inline element_type minCoeff() const { return std::min({x[0], x[1], x[2], x[3]}); }

    // Operators
    inline element_type operator[](int idx) const { return x[idx]; }
    inline element_type& operator[](int idx) { return x[idx]; }
    inline operator const element_type*(void) const { return x; }
    inline Vec4& operator=(const Vec4& v) { x[0]= v.x[0]; x[1]= v.x[1]; x[2]= v.x[2]; x[3]= v.x[3]; return *this; }
    inline bool operator==(const Vec4& v) const { return ((x[0] == v.x[0]) && (x[1] == v.x[1]) && (x[2] == v.x[2]) && (x[3] == v.x[3])); }
    inline Vec4& operator+=(const Vec4& v) { x[0]+= v.x[0]; x[1]+= v.x[1]; x[2]+= v.x[2]; x[3]+= v.x[3]; return *this; }
    inline Vec4& operator-=(const Vec4& v) { x[0]-= v.x[0]; x[1]-= v.x[1]; x[2]-= v.x[2]; x[3]-= v.x[3]; return *this; }
    inline Vec4& operator*=(const element_type f) { x[0]*= f; x[1]*= f; x[2]*= f; x[3]*= f; return *this; }
    inline Vec4& operator/=(const element_type f) { x[0]/= f; x[1]/= f; x[2]/= f; x[3]/= f; return *this; }
    friend Vec4 operator+(Vec4 v, const Vec4& w) noexcept { return v+= w; }
    friend Vec4 operator-(Vec4 v, const Vec4& w) noexcept { return v-= w; }
    friend Vec4 operator*(Vec4 v, const element_type f) noexcept { return v*= f; }
    friend Vec4 operator*(const element_type f, Vec4 v) noexcept { return v*= f; }
    friend Vec4 operator/(Vec4 v, const element_type f) noexcept { return v/= f; }

    inline element_type* array() { return x; }

    private:
    element_type x[4];
  };

}  // namespace Vec
