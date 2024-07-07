#pragma once

// Standard lib
#include <vector>

namespace Field {

  // TODO Remove all nested array functions once nested arrays are swapped for flat arrays everywhere
  // Allocation of field with copy constructor
  template <typename element_type>
  inline std::vector<std::vector<element_type>> AllocNested2(int const iNbA, int const iNbB, element_type const& iVal) {
    return std::vector<std::vector<element_type>>(iNbA, std::vector<element_type>(iNbB, iVal));
  }
  // Allocation of field with copy constructor
  template <typename element_type>
  inline std::vector<std::vector<std::vector<element_type>>> AllocNested3(int const iNbA, int const iNbB, int const iNbC, element_type const& iVal) {
    return std::vector<std::vector<std::vector<element_type>>>(iNbA, std::vector<std::vector<element_type>>(iNbB, std::vector<element_type>(iNbC, iVal)));
  }
  // Allocation of field with copy constructor
  template <typename element_type>
  inline std::vector<std::vector<std::vector<std::vector<element_type>>>> AllocNested4(int const iNbA, int const iNbB, int const iNbC, int const iNbD, element_type const& iVal) {
    return std::vector<std::vector<std::vector<std::vector<element_type>>>>(iNbA, std::vector<std::vector<std::vector<element_type>>>(iNbB, std::vector<std::vector<element_type>>(iNbC, std::vector<element_type>(iNbD, iVal))));
  }
  // Allocation of field with copy constructor
  template <typename element_type>
  inline std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> AllocNested5(int const iNbA, int const iNbB, int const iNbC, int const iNbD, int const iNbE, element_type const& iVal) {
    return std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>>(iNbA, std::vector<std::vector<std::vector<std::vector<element_type>>>>(iNbB, std::vector<std::vector<std::vector<element_type>>>(iNbC, std::vector<std::vector<element_type>>(iNbD, std::vector<element_type>(iNbE, iVal)))));
  }


  // Get dimensions of field
  template <typename element_type>
  inline void GetDim(std::vector<std::vector<element_type>> const& iField, int& oNbA, int& oNbB) {
    oNbA= oNbB= 0;
    oNbA= int(iField.size());
    if (oNbA > 0) {
      oNbB= int(iField[0].size());
    }
  }
  // Get dimensions of field
  template <typename element_type>
  inline void GetDim(std::vector<std::vector<std::vector<element_type>>> const& iField, int& oNbA, int& oNbB, int& oNbC) {
    oNbA= oNbB= oNbC= 0;
    oNbA= int(iField.size());
    if (oNbA > 0) {
      oNbB= int(iField[0].size());
      if (oNbB > 0) {
        oNbC= int(iField[0][0].size());
      }
    }
  }
  // Get dimensions of field
  template <typename element_type>
  inline void GetDim(std::vector<std::vector<std::vector<std::vector<element_type>>>> const& iField, int& oNbA, int& oNbB, int& oNbC, int& oNbD) {
    oNbA= oNbB= oNbC= oNbD= 0;
    oNbA= int(iField.size());
    if (oNbA > 0) {
      oNbB= int(iField[0].size());
      if (oNbB > 0) {
        oNbC= int(iField[0][0].size());
        if (oNbC > 0) {
          oNbD= int(iField[0][0][0].size());
        }
      }
    }
  }
  // Get dimensions of field
  template <typename element_type>
  inline void GetDim(std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> const& iField, int& oNbA, int& oNbB, int& oNbC, int& oNbD, int& oNbE) {
    oNbA= oNbB= oNbC= oNbD= oNbE= 0;
    oNbA= int(iField.size());
    if (oNbA > 0) {
      oNbB= int(iField[0].size());
      if (oNbB > 0) {
        oNbC= int(iField[0][0].size());
        if (oNbC > 0) {
          oNbD= int(iField[0][0][0].size());
          if (oNbD > 0) {
            oNbE= int(iField[0][0][0][0].size());
          }
        }
      }
    }
  }


  // Check if a field has the given dimensions
  template <typename element_type>
  inline bool CheckDim(std::vector<std::vector<element_type>> const& iFieldA, int const iNbA, int const iNbB) {
    if ((int)iFieldA.size() == iNbA) {
      if (iFieldA.empty()) return true;
      if ((int)iFieldA[0].size() == iNbB) {
        return true;
      }
    }
    return false;
  }
  // Check if a field has the given dimensions
  template <typename element_type>
  inline bool CheckDim(std::vector<std::vector<std::vector<element_type>>> const& iFieldA, int const iNbA, int const iNbB, int const iNbC) {
    if ((int)iFieldA.size() == iNbA) {
      if (iFieldA.empty()) return true;
      if ((int)iFieldA[0].size() == iNbB) {
        if (iFieldA[0].empty()) return true;
        if ((int)iFieldA[0][0].size() == iNbC) {
          return true;
        }
      }
    }
    return false;
  }
  // Check if a field has the given dimensions
  template <typename element_type>
  inline bool CheckDim(std::vector<std::vector<std::vector<std::vector<element_type>>>> const& iFieldA, int const iNbA, int const iNbB, int const iNbC, int const iNbD) {
    if ((int)iFieldA.size() == iNbA) {
      if (iFieldA.empty()) return true;
      if ((int)iFieldA[0].size() == iNbB) {
        if (iFieldA[0].empty()) return true;
        if ((int)iFieldA[0][0].size() == iNbC) {
          if (iFieldA[0][0].empty()) return true;
          if ((int)iFieldA[0][0][0].size() == iNbD) {
            return true;
          }
        }
      }
    }
    return false;
  }
  // Check if a field has the given dimensions
  template <typename element_type>
  inline bool CheckDim(std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> const& iFieldA, int const iNbA, int const iNbB, int const iNbC, int const iNbD, int const iNbE) {
    if ((int)iFieldA.size() == iNbA) {
      if (iFieldA.empty()) return true;
      if ((int)iFieldA[0].size() == iNbB) {
        if (iFieldA[0].empty()) return true;
        if ((int)iFieldA[0][0].size() == iNbC) {
          if (iFieldA[0][0].empty()) return true;
          if ((int)iFieldA[0][0][0].size() == iNbD) {
            if (iFieldA[0][0][0].empty()) return true;
            if ((int)iFieldA[0][0][0][0].size() == iNbE) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }


  // Class for storing a field of values
  // All internal data is public to enable optimization at the user's discretion
  // Cannot use element_type = bool due to weird quirk of std::vector
  //
  // Example:
  // Field::Field2<float> fieldA(515, 3, 5.0f);
  // Field::Field2<float> fieldB= fieldA;
  // fieldB= Field::Field2<float>(100, 100, 10.0f);
  //
  // for (int xy= 0; xy < fieldB.nXY; xy++)
  //   fieldB.at(xy)*= 3;
  //
  // for (int x= 0; x < fieldB.nX; x++)
  //   for (int y= 0; y < fieldB.nY; y++)
  //     fieldA.at(x, y)*= 0.5;
  //
  // fieldA.swap(fieldB);
  template <typename element_type>
  class Field2
  {
public:
    int nX, nY;
    int nXY;
    int c0;
    std::vector<element_type> data;

public:
    Field2() {
      nX= nY= 0;
      nXY= 0;
      c0= 0;
    }
    Field2(int const nX, int const nY) {
      if (nX <= 0 || nY <= 0) throw;
      this->nX= nX;
      this->nY= nY;
      this->nXY= nX * nY;
      this->c0= nY;
      this->data= std::vector<element_type>(this->nXY);
    }
    Field2(int const nX, int const nY, element_type const& val) {
      if (nX <= 0 || nY <= 0) throw;
      this->nX= nX;
      this->nY= nY;
      this->nXY= nX * nY;
      this->c0= nY;
      this->data= std::vector<element_type>(this->nXY, val);
    }
    Field2(Field2 const& refField) {
      this->nX= refField.nX;
      this->nY= refField.nY;
      this->nXY= refField.nXY;
      this->c0= refField.c0;
      this->data= refField.data;
    }

    // Swap operator O(1)
    inline void swap(Field2& refField) {
      std::swap(this->nX, refField.nX);
      std::swap(this->nY, refField.nY);
      std::swap(this->nXY, refField.nXY);
      std::swap(this->c0, refField.c0);
      this->data.swap(refField.data);
    }

    // Accessor
    inline element_type& at(int const xy) { return data[xy]; }
    inline const element_type& at(int const xy) const { return data[xy]; }
    inline element_type& at(int const x, int const y) { return data[x * c0 + y]; }
    inline const element_type& at(int const x, int const y) const { return data[x * c0 + y]; }

    // Overload of operator=
    Field2& operator=(Field2 const& refField) {
      if (this == &refField) return *this;
      this->nX= refField.nX;
      this->nY= refField.nY;
      this->nXY= refField.nXY;
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
  // Field::Field3<float> fieldA(515, 3, 45, 5.0f);
  // Field::Field3<float> fieldB= fieldA;
  // fieldB= Field::Field3<float>(100, 100, 100, 10.0f);
  //
  // for (int xyz= 0; xyz < fieldB.nXYZ; xyz++)
  //   fieldB.at(xyz)*= 3;
  //
  // for (int x= 0; x < fieldB.nX; x++)
  //   for (int y= 0; y < fieldB.nY; y++)
  //     for (int z= 0; z < fieldB.nZ; z++)
  //       fieldA.at(x, y, z)*= 0.5;
  //
  // fieldA.swap(fieldB);
  template <typename element_type>
  class Field3
  {
public:
    int nX, nY, nZ;
    int nXYZ;
    int c0, c1;
    std::vector<element_type> data;

public:
    Field3() {
      nX= nY= nZ= 0;
      nXYZ= 0;
      c0= c1= 0;
    }
    Field3(int const nX, int const nY, int const nZ) {
      if (nX <= 0 || nY <= 0 || nZ <= 0) throw;
      this->nX= nX;
      this->nY= nY;
      this->nZ= nZ;
      this->nXYZ= nX * nY * nZ;
      this->c0= nY * nZ;
      this->c1= nZ;
      this->data= std::vector<element_type>(this->nXYZ);
    }
    Field3(int const nX, int const nY, int const nZ, element_type const& val) {
      if (nX <= 0 || nY <= 0 || nZ <= 0) throw;
      this->nX= nX;
      this->nY= nY;
      this->nZ= nZ;
      this->nXYZ= nX * nY * nZ;
      this->c0= nY * nZ;
      this->c1= nZ;
      this->data= std::vector<element_type>(this->nXYZ, val);
    }
    Field3(Field3 const& refField) {
      this->nX= refField.nX;
      this->nY= refField.nY;
      this->nZ= refField.nZ;
      this->nXYZ= refField.nXYZ;
      this->c0= refField.c0;
      this->c1= refField.c1;
      this->data= refField.data;
    }

    // Swap operator O(1)
    inline void swap(Field3& refField) {
      std::swap(this->nX, refField.nX);
      std::swap(this->nY, refField.nY);
      std::swap(this->nZ, refField.nZ);
      std::swap(this->nXYZ, refField.nXYZ);
      std::swap(this->c0, refField.c0);
      std::swap(this->c1, refField.c1);
      this->data.swap(refField.data);
    }

    // Accessor
    inline element_type& at(int const xyz) { return data[xyz]; }
    inline const element_type& at(int const xyz) const { return data[xyz]; }
    inline element_type& at(int const x, int const y, int const z) { return data[x * c0 + y * c1 + z]; }
    inline const element_type& at(int const x, int const y, int const z) const { return data[x * c0 + y * c1 + z]; }

    // Overload of operator=
    Field3& operator=(Field3 const& refField) {
      if (this == &refField) return *this;
      this->nX= refField.nX;
      this->nY= refField.nY;
      this->nZ= refField.nZ;
      this->nXYZ= refField.nXYZ;
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
  // Field::Field4<float> fieldA(2, 515, 3, 45, 5.0f);
  // Field::Field4<float> fieldB= fieldA;
  // fieldB= Field::Field4<float>(100, 100, 100, 2, 10.0f);
  //
  // for (int xyzw= 0; xyzw < fieldB.nXYZW; xyzw++)
  //   fieldB.at(xyzw)*= 3;
  //
  // for (int x= 0; x < fieldB.nX; x++)
  //   for (int y= 0; y < fieldB.nY; y++)
  //     for (int z= 0; z < fieldB.nZ; z++)
  //       for (int w= 0; w < fieldB.nW; w++)
  //         fieldA.at(x, y, z, w)*= 0.5;
  //
  // fieldA.swap(fieldB);
  template <typename element_type>
  class Field4
  {
public:
    int nX, nY, nZ, nW;
    int nXYZW;
    int c0, c1, c2;
    std::vector<element_type> data;

public:
    Field4() {
      nX= nY= nZ= nW= 0;
      nXYZW= 0;
      c0= c1= c2= 0;
    }
    Field4(int const nX, int const nY, int const nZ, int const nW) {
      if (nX <= 0 || nY <= 0 || nZ <= 0 || nW <= 0) throw;
      this->nX= nX;
      this->nY= nY;
      this->nZ= nZ;
      this->nW= nW;
      this->nXYZW= nX * nY * nZ * nW;
      this->c0= nY * nZ * nW;
      this->c1= nZ * nW;
      this->c2= nW;
      this->data= std::vector<element_type>(this->nXYZW);
    }
    Field4(int const nX, int const nY, int const nZ, int const nW, element_type const& val) {
      if (nX <= 0 || nY <= 0 || nZ <= 0 || nW <= 0) throw;
      this->nX= nX;
      this->nY= nY;
      this->nZ= nZ;
      this->nW= nW;
      this->nXYZW= nX * nY * nZ * nW;
      this->c0= nY * nZ * nW;
      this->c1= nZ * nW;
      this->c2= nW;
      this->data= std::vector<element_type>(this->nXYZW, val);
    }
    Field4(Field4 const& refField) {
      this->nX= refField.nX;
      this->nY= refField.nY;
      this->nZ= refField.nZ;
      this->nW= refField.nW;
      this->nXYZW= refField.nXYZW;
      this->c0= refField.c0;
      this->c1= refField.c1;
      this->c2= refField.c2;
      this->data= refField.data;
    }

    // Swap operator O(1)
    inline void swap(Field4& refField) {
      std::swap(this->nX, refField.nX);
      std::swap(this->nY, refField.nY);
      std::swap(this->nZ, refField.nZ);
      std::swap(this->nW, refField.nW);
      std::swap(this->nXYZW, refField.nXYZW);
      std::swap(this->c0, refField.c0);
      std::swap(this->c1, refField.c1);
      std::swap(this->c2, refField.c2);
      this->data.swap(refField.data);
    }

    // Accessor
    inline element_type& at(int const xyzw) { return data[xyzw]; }
    inline const element_type& at(int const xyzw) const { return data[xyzw]; }
    inline element_type& at(int const x, int const y, int const z, int const w) { return data[x * c0 + y * c1 + z * c2 + w]; }
    inline const element_type& at(int const x, int const y, int const z, int const w) const { return data[x * c0 + y * c1 + z * c2 + w]; }

    // Overload of operator=
    Field4& operator=(Field4 const& refField) {
      if (this == &refField) return *this;
      this->nX= refField.nX;
      this->nY= refField.nY;
      this->nZ= refField.nZ;
      this->nW= refField.nW;
      this->nXYZW= refField.nXYZW;
      this->c0= refField.c0;
      this->c1= refField.c1;
      this->c2= refField.c2;
      this->data= refField.data;
      return *this;
    }
  };
}  // namespace Field
