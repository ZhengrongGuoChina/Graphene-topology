// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGYVECTOR3HPP
#define GRAPOLOGYVECTOR3HPP

namespace grapology {

using DataTypeInt = int;
using DataTypeFloat = float;

/**
 * @brief A template class representing a 3D vector or point in space
 *
 * This class provides complete 3D vector operations including basic arithmetic,
 * dot product, cross product, normalization, rotation, and comparison
 * operations. As a template class, it can accept any arithmetic type as the
 * vector element type (e.g., int, float, double).
 *
 * @tparam TI The type of the vector elements, must satisfy the VecEleType
 * concept (arithmetic type)
 */

template <typename T>
concept VecEleType =
    std::is_arithmetic_v<T>;  // Vector Element Type must be arithmetic type

template <VecEleType TI>
struct vector3 final {
  using value_type = TI;

  value_type x, y, z;

  constexpr explicit vector3() : x(0), y(0), z(0) {}

  constexpr vector3(const value_type ixd,
                    const value_type iyd,
                    const value_type izd)
      : x(ixd), y(iyd), z(izd) {}

  template <VecEleType TO>
  constexpr vector3(const TO ixd,
                    const TO iyd,
                    const TO izd)
      : x(static_cast<value_type>(ixd)), y(static_cast<value_type>(iyd)), z(static_cast<value_type>(izd)) {}

  constexpr vector3(const vector3& ru) : x(ru.x), y(ru.y), z(ru.z) {}

  template <VecEleType TO>
  constexpr vector3(const vector3<TO>& ru)
      : x(static_cast<value_type>(ru.x)),
        y(static_cast<value_type>(ru.y)),
        z(static_cast<value_type>(ru.z)) {}

  template <VecEleType TO>
  constexpr vector3(std::initializer_list<TO> args):
   x(static_cast<value_type>(*args.begin())),
   y(static_cast<value_type>(*(args.begin() + 1))),
   z(static_cast<value_type>(*(args.begin() + 2))) {
    //assert(args.size() == 3, "vector3 requires 3 elements");
  }
  
  constexpr auto&& operator[](const size_t index) {
    switch (index) {
      case 0:
        return x;
      case 1:
        return y;
      default:
        return z;
    }
  }

  constexpr auto operator[](const size_t index) const -> const value_type {
    switch (index) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default: {
        //assert(index < 3, "vector3 index overflow");
        return x;  // less warnings
      }
    }
  }

  template <VecEleType TO>
  constexpr inline auto operator=(const vector3<TO>& ru) noexcept -> vector3& {
    x = static_cast<value_type>(ru.x);
    y = static_cast<value_type>(ru.y);
    z = static_cast<value_type>(ru.z);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator+=(const vector3<TO>& ru) noexcept -> vector3& {
    x += static_cast<value_type>(ru.x);
    y += static_cast<value_type>(ru.y);
    z += static_cast<value_type>(ru.z);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator+=(const TO ru) noexcept -> vector3& {
    x += static_cast<value_type>(ru);
    y += static_cast<value_type>(ru);
    z += static_cast<value_type>(ru);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator-=(const vector3<TO>& ru) noexcept -> vector3& {
    x -= static_cast<value_type>(ru.x);
    y -= static_cast<value_type>(ru.y);
    z -= static_cast<value_type>(ru.z);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator-=(const TO ru) noexcept -> vector3& {
    x -= static_cast<value_type>(ru);
    y -= static_cast<value_type>(ru);
    z -= static_cast<value_type>(ru);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator*=(const TO ru) noexcept -> vector3& {
    x *= static_cast<value_type>(ru);
    y *= static_cast<value_type>(ru);
    z *= static_cast<value_type>(ru);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator/=(const TO ru) noexcept -> vector3& {
    x /= static_cast<value_type>(ru);
    y /= static_cast<value_type>(ru);
    z /= static_cast<value_type>(ru);
    return *this;
  }

  template <VecEleType TO>
  constexpr inline auto operator+(const vector3<TO>& ru) const noexcept {
    return vector3<decltype(TO() + value_type())>(x + ru.x, y + ru.y, z + ru.z);
  }

  template <VecEleType TO>
  constexpr inline auto operator+(const TO ru) const noexcept {
    return vector3<decltype(TO() + value_type())>(x + ru, y + ru, z + ru);
  }

  template <VecEleType TO>
  constexpr inline auto operator-(const vector3<TO>& ru) const noexcept {
    return vector3<decltype(TO() - value_type())>(x - ru.x, y - ru.y, z - ru.z);
  }

  template <VecEleType TO>
  constexpr inline auto operator-(const TO ru) const noexcept {
    return vector3<decltype(TO() - value_type())>(x - ru, y - ru, z - ru);
  }

  constexpr inline auto operator-(int) const noexcept -> vector3 {
    return {-x, -y, -z};
  }

  template <VecEleType TO>
  constexpr inline auto operator*(const vector3<TO>& ru) const noexcept {
    return x * ru.x + y * ru.y + z * ru.z;
  }

  template <VecEleType TO>
  constexpr inline auto operator*(const TO ru) const noexcept {
    return vector3<decltype(TO() * value_type())>(x * ru, y * ru, z * ru);
  }

  template <VecEleType TO>
  constexpr inline auto operator/(const TO ru) const noexcept {
    return vector3<decltype(TO() * value_type())>(x / ru, y / ru, z / ru);
  }

  template <VecEleType TO>
  constexpr inline auto operator^(const vector3<TO>& ru) const noexcept {
    return vector3<decltype(TO() * value_type())>(
        y * ru.z - z * ru.y, z * ru.x - x * ru.z, x * ru.y - y * ru.x);
  }

  template <VecEleType TO>
  constexpr inline bool operator==(const vector3<TO>& ru) const {
    if constexpr (std::is_integral_v<value_type> && std::is_integral_v<TO>) {
      return x == ru.x && y == ru.y && z == ru.z;
    } else {
      //error(np_debug_info,
      static_assert(false,"Can not compare equal or unequal between floating point vector3s");
      return false;  // less warnings
    }
  }

  template <VecEleType TO>
  constexpr inline bool operator!=(const vector3<TO>& ru) const {
    return not(*this == ru);
  }

  constexpr inline auto pow(size_t ex) const noexcept {
    return vector3(std::pow(x, ex), std::pow(y, ex), std::pow(z, ex));
  }

  constexpr inline auto normal() const noexcept -> vector3<DataTypeFloat> {
    if constexpr (std::is_integral_v<value_type>) {
      DataTypeFloat _nl =
          ::sqrt(static_cast<DataTypeFloat>(x * x + y * y + z * z));
      if (_nl == 0)
        return vector3<DataTypeFloat>(0, 0, 0);
      else {
        const DataTypeFloat _nlq = static_cast<DataTypeFloat>(1.0 / _nl);
        return vector3<DataTypeFloat>(x * _nlq, y * _nlq, z * _nlq);
      }
    } else {
      const value_type _nl = ::sqrt(x * x + y * y + z * z);
      if (_nl == 0)
        return vector3<DataTypeFloat>(0, 0, 0);
      else {
        const DataTypeFloat _nlq = static_cast<DataTypeFloat>(1.0 / _nl);
        return vector3<DataTypeFloat>(x * _nlq, y * _nlq, z * _nlq);
      }
    }
  }

  constexpr inline void normalize() noexcept {
    if constexpr (std::is_integral_v<value_type>) {
      static_assert(false,"Can not normalize an integral vector3");
    } else {
      const value_type _nl = ::sqrt(x * x + y * y + z * z);
      if (_nl != 0) {
        const value_type _nlq = static_cast<value_type>(1.0 / _nl);
        x *= _nlq;
        y *= _nlq;
        z *= _nlq;
      }
    }
  }

  constexpr inline auto length() const noexcept -> DataTypeFloat {
    return ::sqrt(static_cast<DataTypeFloat>(x * x + y * y + z * z));
  }

  constexpr inline auto norm() const noexcept { return x * x + y * y + z * z; }

  constexpr inline void rotate(const vector3& axis,
                               const DataTypeFloat angle) noexcept {
    vector3 ax = axis.normal();

    DataTypeFloat matrix[3][3];

    matrix[0][0] = ::cos(angle) + (ax.x * ax.x) * (1 - ::cos(angle));
    matrix[0][1] = (ax.x * ax.y) * (1 - ::cos(angle)) - ax.z * ::sin(angle);
    matrix[0][2] = (ax.x * ax.z) * (1 - ::cos(angle)) + ax.y * ::sin(angle);

    matrix[1][0] = (ax.x * ax.y) * (1 - ::cos(angle)) + ax.z * ::sin(angle);
    matrix[1][1] = ::cos(angle) + (ax.y * ax.y) * (1 - ::cos(angle));
    matrix[1][2] = (ax.y * ax.z) * (1 - ::cos(angle)) - ax.x * ::sin(angle);

    matrix[2][0] = (ax.x * ax.z) * (1 - ::cos(angle)) - ax.y * ::sin(angle);
    matrix[2][1] = (ax.y * ax.z) * (1 - ::cos(angle)) + ax.x * ::sin(angle);
    matrix[2][2] = ::cos(angle) + (ax.z * ax.z) * (1 - ::cos(angle));

    value_type _x = static_cast<value_type>(
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z);
    value_type _y = static_cast<value_type>(
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z);
    value_type _z = static_cast<value_type>(
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z);

    x = _x;
    y = _y;
    z = _z;
  }

  template <VecEleType TO>
  constexpr inline auto angle(const vector3<TO>& axis) const noexcept -> DataTypeFloat {
    return ::acos(static_cast<DataTypeFloat>(
        (x * axis.x + y * axis.y + z * axis.z) /
        ::sqrt((x * x + y * y + z * z) *
               (axis.x * axis.x + axis.y * axis.y + axis.z * axis.z))));
  }

};  // end of class vector3

template <VecEleType TI, VecEleType TO>
constexpr inline auto operator+(const TO r, const vector3<TI>& ru) noexcept {
  return vector3<decltype(TO() + TI())>(r + ru.x, r + ru.y, r + ru.z);
}

template <VecEleType TI, VecEleType TO>
constexpr inline auto operator-(const TO r, const vector3<TI>& ru) noexcept {
  return vector3<decltype(TO() - TI())>(r - ru.x, r - ru.y, r - ru.z);
}

template <VecEleType TI, VecEleType TO>
constexpr inline auto operator*(const TO r, const vector3<TI>& ru) noexcept {
  return vector3<decltype(TO() * TI())>(r * ru.x, r * ru.y, r * ru.z);
}

// only use this to create reciprocal vector3
template <VecEleType TI, VecEleType TO>
constexpr inline auto operator/(const TO r, const vector3<TI>& ru) noexcept {
  return vector3<decltype(TO() / TI())>(r / ru.x, r / ru.y, r / ru.z);
}


using DataTypeVector3D = vector3<DataTypeFloat>;

}  // namespace grapology

#endif
