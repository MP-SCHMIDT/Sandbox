#pragma once

// Standard lib
#include <vector>


// Class for storing a field of values
// All internal data is public to enable optimization at the user's discretion
// Cannot use element_type = bool due to weird quirk of std::vector
//
// Example:
// Field2D<float> fieldA(515, 3, 5.0f);
// Field2D<float> fieldB= fieldA;
// fieldB= Field2D<float>(100, 100, 10.0f);
//
// for (int k= 0; k < fieldB.nK; k++)
//   fieldB(k)*= 3;
//
// for (int x= 0; x < fieldB.nX; x++)
//   for (int y= 0; y < fieldB.nY; y++)
//     fieldA(x, y)*= 0.5;
//
// fieldA.swap(fieldB);
template <typename element_type>
class Field2D
{
  public:
  int nX, nY;
  int nK;
  int c0;
  std::vector<element_type> data;

  public:
  // Constructors
  Field2D() {
    nX= nY= 0;
    nK= 0;
    c0= 0;
  }

  Field2D(int const nX, int const nY) {
    if (nX <= 0 || nY <= 0) throw;
    this->nX= nX;
    this->nY= nY;
    this->nK= nX * nY;
    this->c0= nY;
    this->data= std::vector<element_type>(this->nK);
  }

  Field2D(int const nX, int const nY, element_type const& val) {
    if (nX <= 0 || nY <= 0) throw;
    this->nX= nX;
    this->nY= nY;
    this->nK= nX * nY;
    this->c0= nY;
    this->data= std::vector<element_type>(this->nK, val);
  }

  Field2D(Field2D const& refField) {
    this->nX= refField.nX;
    this->nY= refField.nY;
    this->nK= refField.nK;
    this->c0= refField.c0;
    this->data= refField.data;
  }

  // Swap operator
  inline void swap(Field2D& refField) {
    std::swap(this->nX, refField.nX);
    std::swap(this->nY, refField.nY);
    std::swap(this->nK, refField.nK);
    std::swap(this->c0, refField.c0);
    this->data.swap(refField.data);
  }

  // Overloads of operator() to access data
  inline element_type& operator()(int const k) { return data[k]; }
  inline const element_type& operator()(int const k) const { return data[k]; }
  inline element_type& operator()(int const x, int const y) { return data[x * c0 + y]; }
  inline const element_type& operator()(int const x, int const y) const { return data[x * c0 + y]; }

  // Overload of operator=
  Field2D& operator=(Field2D const& refField) {
    if (this == &refField) return *this;
    this->nX= refField.nX;
    this->nY= refField.nY;
    this->nK= refField.nK;
    this->c0= refField.c0;
    this->data= refField.data;
    return *this;
  }
};


// Class for storing a field of values
// All internal data is public to enable optimization at the user's discretion
// Cannot use element_type = bool due to weird quirk of std::vector
//
// Example:
// Field3D<float> fieldA(515, 3, 45, 5.0f);
// Field3D<float> fieldB= fieldA;
// fieldB= Field3D<float>(100, 100, 100, 10.0f);
//
// for (int k= 0; k < fieldB.nK; k++)
//   fieldB(k)*= 3;
//
// for (int x= 0; x < fieldB.nX; x++)
//   for (int y= 0; y < fieldB.nY; y++)
//     for (int z= 0; z < fieldB.nZ; z++)
//       fieldA(x, y, z)*= 0.5;
//
// fieldA.swap(fieldB);
template <typename element_type>
class Field3D
{
  public:
  int nX, nY, nZ;
  int nK;
  int c0, c1;
  std::vector<element_type> data;

  public:
  // Constructors
  Field3D() {
    nX= nY= nZ= 0;
    nK= 0;
    c0= c1= 0;
  }

  Field3D(int const nX, int const nY, int const nZ) {
    if (nX <= 0 || nY <= 0 || nZ <= 0) throw;
    this->nX= nX;
    this->nY= nY;
    this->nZ= nZ;
    this->nK= nX * nY * nZ;
    this->c0= nY * nZ;
    this->c1= nZ;
    this->data= std::vector<element_type>(this->nK);
  }

  Field3D(int const nX, int const nY, int const nZ, element_type const& val) {
    if (nX <= 0 || nY <= 0 || nZ <= 0) throw;
    this->nX= nX;
    this->nY= nY;
    this->nZ= nZ;
    this->nK= nX * nY * nZ;
    this->c0= nY * nZ;
    this->c1= nZ;
    this->data= std::vector<element_type>(this->nK, val);
  }

  Field3D(Field3D const& refField) {
    this->nX= refField.nX;
    this->nY= refField.nY;
    this->nZ= refField.nZ;
    this->nK= refField.nK;
    this->c0= refField.c0;
    this->c1= refField.c1;
    this->data= refField.data;
  }

  // Swap operator
  inline void swap(Field3D& refField) {
    std::swap(this->nX, refField.nX);
    std::swap(this->nY, refField.nY);
    std::swap(this->nZ, refField.nZ);
    std::swap(this->nK, refField.nK);
    std::swap(this->c0, refField.c0);
    std::swap(this->c1, refField.c1);
    this->data.swap(refField.data);
  }

  // Overloads of operator() to access data
  inline element_type& operator()(int const k) { return data[k]; }
  inline const element_type& operator()(int const k) const { return data[k]; }
  inline element_type& operator()(int const x, int const y, int const z) { return data[x * c0 + y * c1 + z]; }
  inline const element_type& operator()(int const x, int const y, int const z) const { return data[x * c0 + y * c1 + z]; }

  // Overload of operator=
  Field3D& operator=(Field3D const& refField) {
    if (this == &refField) return *this;
    this->nX= refField.nX;
    this->nY= refField.nY;
    this->nZ= refField.nZ;
    this->nK= refField.nK;
    this->c0= refField.c0;
    this->c1= refField.c1;
    this->data= refField.data;
    return *this;
  }
};


// Class for storing a field of values
// All internal data is public to enable optimization at the user's discretion
// Cannot use element_type = bool due to weird quirk of std::vector
//
// Example:
// Field4D<float> fieldA(2, 515, 3, 45, 5.0f);
// Field4D<float> fieldB= fieldA;
// fieldB= Field4D<float>(100, 100, 100, 2, 10.0f);
//
// for (int k= 0; k < fieldB.nK; k++)
//   fieldB(k)*= 3;
//
// for (int x= 0; x < fieldB.nX; x++)
//   for (int y= 0; y < fieldB.nY; y++)
//     for (int z= 0; z < fieldB.nZ; z++)
//       for (int w= 0; w < fieldB.nW; w++)
//         fieldA(x, y, z, w)*= 0.5;
//
// fieldA.swap(fieldB);
template <typename element_type>
class Field4D
{
  public:
  int nX, nY, nZ, nW;
  int nK;
  int c0, c1, c2;
  std::vector<element_type> data;

  public:
  // Constructors
  Field4D() {
    nX= nY= nZ= nW= 0;
    nK= 0;
    c0= c1= c2= 0;
  }

  Field4D(int const nX, int const nY, int const nZ, int const nW) {
    if (nX <= 0 || nY <= 0 || nZ <= 0 || nW <= 0) throw;
    this->nX= nX;
    this->nY= nY;
    this->nZ= nZ;
    this->nW= nW;
    this->nK= nX * nY * nZ * nW;
    this->c0= nY * nZ * nW;
    this->c1= nZ * nW;
    this->c2= nW;
    this->data= std::vector<element_type>(this->nK);
  }

  Field4D(int const nX, int const nY, int const nZ, int const nW, element_type const& val) {
    if (nX <= 0 || nY <= 0 || nZ <= 0 || nW <= 0) throw;
    this->nX= nX;
    this->nY= nY;
    this->nZ= nZ;
    this->nW= nW;
    this->nK= nX * nY * nZ * nW;
    this->c0= nY * nZ * nW;
    this->c1= nZ * nW;
    this->c2= nW;
    this->data= std::vector<element_type>(this->nK, val);
  }

  Field4D(Field4D const& refField) {
    this->nX= refField.nX;
    this->nY= refField.nY;
    this->nZ= refField.nZ;
    this->nW= refField.nW;
    this->nK= refField.nK;
    this->c0= refField.c0;
    this->c1= refField.c1;
    this->c2= refField.c2;
    this->data= refField.data;
  }

  // Swap operator
  inline void swap(Field4D& refField) {
    std::swap(this->nX, refField.nX);
    std::swap(this->nY, refField.nY);
    std::swap(this->nZ, refField.nZ);
    std::swap(this->nW, refField.nW);
    std::swap(this->nK, refField.nK);
    std::swap(this->c0, refField.c0);
    std::swap(this->c1, refField.c1);
    std::swap(this->c2, refField.c2);
    this->data.swap(refField.data);
  }

  // Overloads of operator() to access data
  inline element_type& operator()(int const k) { return data[k]; }
  inline const element_type& operator()(int const k) const { return data[k]; }
  inline element_type& operator()(int const x, int const y, int const z, int const w) { return data[x * c0 + y * c1 + z * c2 + w]; }
  inline const element_type& operator()(int const x, int const y, int const z, int const w) const { return data[x * c0 + y * c1 + z * c2 + w]; }

  // Overload of operator=
  Field4D& operator=(Field4D const& refField) {
    if (this == &refField) return *this;
    this->nX= refField.nX;
    this->nY= refField.nY;
    this->nZ= refField.nZ;
    this->nW= refField.nW;
    this->nK= refField.nK;
    this->c0= refField.c0;
    this->c1= refField.c1;
    this->c2= refField.c2;
    this->data= refField.data;
    return *this;
  }
};
