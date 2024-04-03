#pragma once

// Standard lib
#include <vector>


namespace Field {

  // Allocation of fields
  template <typename element_type>
  inline std::vector<std::vector<element_type>> AllocField2D(int const iNbA, int const iNbB, element_type const& val) {
    return std::vector<std::vector<element_type>>(iNbA, std::vector<element_type>(iNbB, val));
  }
  template <typename element_type>
  inline std::vector<std::vector<std::vector<element_type>>> AllocField3D(int const iNbA, int const iNbB, int const iNbC, element_type const& val) {
    return std::vector<std::vector<std::vector<element_type>>>(iNbA, std::vector<std::vector<element_type>>(iNbB, std::vector<element_type>(iNbC, val)));
  }
  template <typename element_type>
  inline std::vector<std::vector<std::vector<std::vector<element_type>>>> AllocField4D(int const iNbA, int const iNbB, int const iNbC, int const iNbD, element_type const& val) {
    return std::vector<std::vector<std::vector<std::vector<element_type>>>>(iNbA, std::vector<std::vector<std::vector<element_type>>>(iNbB, std::vector<std::vector<element_type>>(iNbC, std::vector<element_type>(iNbD, val))));
  }
  template <typename element_type>
  inline std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> AllocField5D(int const iNbA, int const iNbB, int const iNbC, int const iNbD, int const iNbE, element_type const& val) {
    return std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>>(iNbA, std::vector<std::vector<std::vector<std::vector<element_type>>>>(iNbB, std::vector<std::vector<std::vector<element_type>>>(iNbC, std::vector<std::vector<element_type>>(iNbD, std::vector<element_type>(iNbE, val)))));
  }

  // Get dimensions of fields
  template <typename element_type>
  inline void GetFieldDimensions(std::vector<std::vector<element_type>> const& iField, int& oNbA, int& oNbB) {
    oNbA= oNbB= 0;
    oNbA= int(iField.size());
    if (oNbA > 0) {
      oNbB= int(iField[0].size());
    }
  }
  template <typename element_type>
  inline void GetFieldDimensions(std::vector<std::vector<std::vector<element_type>>> const& iField, int& oNbA, int& oNbB, int& oNbC) {
    oNbA= oNbB= oNbC= 0;
    oNbA= int(iField.size());
    if (oNbA > 0) {
      oNbB= int(iField[0].size());
      if (oNbB > 0) {
        oNbC= int(iField[0][0].size());
      }
    }
  }
  template <typename element_type>
  inline void GetFieldDimensions(std::vector<std::vector<std::vector<std::vector<element_type>>>> const& iField, int& oNbA, int& oNbB, int& oNbC, int& oNbD) {
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
  template <typename element_type>
  inline void GetFieldDimensions(std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> const& iField, int& oNbA, int& oNbB, int& oNbC, int& oNbD, int& oNbE) {
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

  // Check if two fields share the same dimensions
  template <typename element_type>
  inline bool CheckSameDimensions(std::vector<std::vector<element_type>> const& iFieldA,
                                  std::vector<std::vector<element_type>> const& iFieldB) {
    if (iFieldA.size() == iFieldB.size()) {
      if (iFieldA.empty()) return true;
      if (iFieldA[0].size() == iFieldB[0].size()) {
        return true;
      }
    }
    return false;
  }
  template <typename element_type>
  inline bool CheckSameDimensions(std::vector<std::vector<std::vector<element_type>>> const& iFieldA,
                                  std::vector<std::vector<std::vector<element_type>>> const& iFieldB) {
    if (iFieldA.size() == iFieldB.size()) {
      if (iFieldA.empty()) return true;
      if (iFieldA[0].size() == iFieldB[0].size()) {
        if (iFieldA[0].empty()) return true;
        if (iFieldA[0][0].size() == iFieldB[0][0].size()) {
          return true;
        }
      }
    }
    return false;
  }
  template <typename element_type>
  inline bool CheckSameDimensions(std::vector<std::vector<std::vector<std::vector<element_type>>>> const& iFieldA,
                                  std::vector<std::vector<std::vector<std::vector<element_type>>>> const& iFieldB) {
    if (iFieldA.size() == iFieldB.size()) {
      if (iFieldA.empty()) return true;
      if (iFieldA[0].size() == iFieldB[0].size()) {
        if (iFieldA[0].empty()) return true;
        if (iFieldA[0][0].size() == iFieldB[0][0].size()) {
          if (iFieldA[0][0].empty()) return true;
          if (iFieldA[0][0][0].size() == iFieldB[0][0][0].size()) {
            return true;
          }
        }
      }
    }
    return false;
  }
  template <typename element_type>
  inline bool CheckSameDimensions(std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> const& iFieldA,
                                  std::vector<std::vector<std::vector<std::vector<std::vector<element_type>>>>> const& iFieldB) {
    if (iFieldA.size() == iFieldB.size()) {
      if (iFieldA.empty()) return true;
      if (iFieldA[0].size() == iFieldB[0].size()) {
        if (iFieldA[0].empty()) return true;
        if (iFieldA[0][0].size() == iFieldB[0][0].size()) {
          if (iFieldA[0][0].empty()) return true;
          if (iFieldA[0][0][0].size() == iFieldB[0][0][0].size()) {
            if (iFieldA[0][0][0].empty()) return true;
            if (iFieldA[0][0][0][0].size() == iFieldB[0][0][0][0].size()) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }
}  // namespace Field


/*
// Standard lib
#include <vector>

namespace Math {

  template <typename element_type>
  class Field1D
  {
public:
    int nbX;
    int nbK;
    std::vector<element_type> data;

public:
    // Constructors
    Field1D() {
      nbX= 0;
      nbK= 0;
    }

    Field1D(int const nbX) {
      if (nbX < 0) throw;
      this->nbX= nbX;
      this->nbK= nbX;
      data= std::vector<element_type>(nbX);
    }

    Field1D(int const nbX, int const nbY, element_type const& val) {
      if (nbX < 0) throw;
      this->nbX= nbX;
      this->nbK= nbX;
      data= std::vector<element_type>(nbX, val);
    }

    Field1D(Field1D const& refField) {
      this->nbX= refField.nbX;
      this->nbK= refField.nbX;
      data= refField;
    }

    // Overloads of operator() to access data
    inline element_type& operator()(int const idxX) {
      return data[idxX];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxX) const {
      return data[idxX];  // Cannot use element_type = bool due to weird quirk of std::vector
    }

    // Overload of operator=
    Field1D& operator=(Field1D const& refField) {
      if (this == &refField)
        return *this;
      this->nbX= refField.nbX;
      this->nbK= refField.nbK;
      data= refField.data;
      return *this;
    }
  };


  template <typename element_type>
  class Field2D
  {
public:
    int nbX, nbY;
    int nbK;
    std::vector<std::vector<element_type>> data;

public:
    // Constructors
    Field2D() {
      nbX= 0;
      nbY= 0;
      nbK= 0;
    }

    Field2D(int const nbX, int const nbY) {
      if (nbX < 0 || nbY < 0) throw;
      this->nbX= nbX;
      this->nbY= nbY;
      this->nbK= nbX * nbY;
      data= std::vector<std::vector<element_type>>(nbX, std::vector<element_type>(nbY));
    }

    Field2D(int const nbX, int const nbY, element_type const& val) {
      if (nbX < 0 || nbY < 0) throw;
      this->nbX= nbX;
      this->nbY= nbY;
      this->nbK= nbX * nbY;
      data= std::vector<std::vector<element_type>>(nbX, std::vector<element_type>(nbY, val));
    }

    Field2D(Field2D const& refField) {
      this->nbX= refField.nbX;
      this->nbY= refField.nbY;
      this->nbK= refField.nbX * refField.nbY;
      data= refField;
    }

    // Overloads of operator() to access data
    inline element_type& operator()(int const idxK) {
      return data[idxK / nbY][idxK % nbY];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxK) const {
      return data[idxK / nbY][idxK % nbY];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline element_type& operator()(int const idxX, int const idxY) {
      return data[idxX][idxY];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxX, int const idxY) const {
      return data[idxX][idxY];  // Cannot use element_type = bool due to weird quirk of std::vector
    }

    // Overload of operator=
    Field2D& operator=(Field2D const& refField) {
      if (this == &refField)
        return *this;
      this->nbX= refField.nbX;
      this->nbY= refField.nbY;
      this->nbK= refField.nbK;
      data= refField.data;
      return *this;
    }
  };


  template <typename element_type>
  class Field3D
  {
public:
    int nbX, nbY, nbZ;
    int nbK;
    std::vector<std::vector<std::vector<element_type>>> data;

public:
    // Constructors
    Field3D() {
      nbX= 0;
      nbY= 0;
      nbZ= 0;
      nbK= 0;
    }

    Field3D(int const nbX, int const nbY, int const nbZ) {
      if (nbX < 0 || nbY < 0 || nbZ < 0) throw;
      this->nbX= nbX;
      this->nbY= nbY;
      this->nbZ= nbZ;
      this->nbK= nbX * nbY * nbZ;
      data= std::vector<std::vector<std::vector<element_type>>>(nbX, std::vector<std::vector<element_type>>(nbY, std::vector<element_type>(nbZ)));
    }

    Field3D(int const nbX, int const nbY, int const nbZ, element_type const& val) {
      if (nbX < 0 || nbY < 0 || nbZ < 0) throw;
      this->nbX= nbX;
      this->nbY= nbY;
      this->nbZ= nbZ;
      this->nbK= nbX * nbY * nbZ;
      data= std::vector<std::vector<std::vector<element_type>>>(nbX, std::vector<std::vector<element_type>>(nbY, std::vector<element_type>(nbZ, val)));
    }

    Field3D(Field3D const& refField) {
      this->nbX= refField.nbX;
      this->nbY= refField.nbY;
      this->nbZ= refField.nbZ;
      this->nbK= refField.nbX * refField.nbY * refField.nbZ;
      data= refField;
    }

    // Overloads of operator() to access data
    inline element_type& operator()(int const idxK) {
      return data[idxK / (nbY * nbZ)][(idxK % nbY) / nbZ][idxK % nbZ];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxK) const {
      return data[idxK / (nbY * nbZ)][(idxK % nbY) / nbZ][idxK % nbZ];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline element_type& operator()(int const idxX, int const idxY, int const idxZ) {
      return data[idxX][idxY][idxZ];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxX, int const idxY, int const idxZ) const {
      return data[idxX][idxY][idxZ];  // Cannot use element_type = bool due to weird quirk of std::vector
    }

    // Overload of operator=
    Field3D& operator=(Field3D const& refField) {
      if (this == &refField)
        return *this;
      this->nbX= refField.nbX;
      this->nbY= refField.nbY;
      this->nbZ= refField.nbZ;
      this->nbK= refField.nbK;
      data= refField.data;
      return *this;
    }
  };


  template <typename element_type>
  class Field4D
  {
public:
    int nbX, nbY, nbZ, nbW;
    int nbK;
    std::vector<std::vector<std::vector<std::vector<element_type>>>> data;

public:
    // Constructors
    Field4D() {
      nbX= 0;
      nbY= 0;
      nbZ= 0;
      nbW= 0;
      nbK= 0;
    }

    Field4D(int const nbX, int const nbY, int const nbZ, int const nbW) {
      if (nbX < 0 || nbY < 0 || nbZ < 0 || nbW < 0) throw;
      this->nbX= nbX;
      this->nbY= nbY;
      this->nbZ= nbZ;
      this->nbW= nbW;
      this->nbK= nbX * nbY * nbZ * nbW;
      data= std::vector<std::vector<std::vector<std::vector<element_type>>>>(nbX, std::vector<std::vector<std::vector<element_type>>>(nbY, std::vector<std::vector<element_type>>(nbZ, std::vector<element_type>(nbW))));
    }

    Field4D(int const nbX, int const nbY, int const nbZ, int const nbW, element_type const& val) {
      if (nbX < 0 || nbY < 0 || nbZ < 0 || nbW < 0) throw;
      this->nbX= nbX;
      this->nbY= nbY;
      this->nbZ= nbZ;
      this->nbW= nbW;
      this->nbK= nbX * nbY * nbZ * nbW;
      data= std::vector<std::vector<std::vector<std::vector<element_type>>>>(nbX, std::vector<std::vector<std::vector<element_type>>>(nbY, std::vector<std::vector<element_type>>(nbZ, std::vector<element_type>(nbW, val))));
    }

    Field4D(Field4D const& refField) {
      this->nbX= refField.nbX;
      this->nbY= refField.nbY;
      this->nbZ= refField.nbZ;
      this->nbW= refField.nbW;
      this->nbK= refField.nbX * refField.nbY * refField.nbZ * refField.nbW;
      data= refField;
    }

    // Overloads of operator() to access data
    inline element_type& operator()(int const idxK) {
      return data[idxK / (nbY * nbZ * nbW)][(idxK % nbY) / (nbZ * nbW)][(idxK % nbZ) / nbW][(idxK % nbW)];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxK) const {
      return data[idxK / (nbY * nbZ * nbW)][(idxK % nbY) / (nbZ * nbW)][(idxK % nbZ) / nbW][(idxK % nbW)];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline element_type& operator()(int const idxX, int const idxY, int const idxZ, int const idxW) {
      return data[idxX][idxY][idxZ][idxW];  // Cannot use element_type = bool due to weird quirk of std::vector
    }
    inline const element_type& operator()(int const idxX, int const idxY, int const idxZ, int const idxW) const {
      return data[idxX][idxY][idxZ][idxW];  // Cannot use element_type = bool due to weird quirk of std::vector
    }

    // Overload of operator=
    Field4D& operator=(Field4D const& refField) {
      if (this == &refField)
        return *this;
      this->nbX= refField.nbX;
      this->nbY= refField.nbY;
      this->nbZ= refField.nbZ;
      this->nbW= refField.nbW;
      this->nbK= refField.nbK;
      data= refField.data;
      return *this;
    }
  };

}  // namespace Math
*/


/*
#include <vector>

template<typename Impl, std::size_t dim>
struct dim_wrapper {
    using type_t = std::vector<typename dim_wrapper<Impl, dim - 1>::type_t>;
};

template<typename Impl>
struct dim_wrapper<Impl, 0> {
    using type_t = Impl;
};

template<typename Impl, std::size_t dim>
using vec_D = typename dim_wrapper<Impl, dim>::type_t;

vec_D<float, 3> myVec;


// using Vec5i= vec_D<int, 5>;
// template<size_t dim> using Vecu= Vec_D<unsigned, dim>;
*/
