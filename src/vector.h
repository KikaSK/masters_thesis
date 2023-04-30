#pragma once

#include <ginac/ginac.h>

#include <iostream>

#include "assertm.h"
#include "constants.h"
#include "point.h"

using GiNaC::numeric;
using std::endl;

class Point;

class Vector {
 private:
  numeric _x, _y, _z;

 public:
  Vector(numeric x, numeric y, numeric z);
  Vector(Point A, Point B);
  Vector(const Vector &v) = default;
  Vector() = delete;

  numeric x() const;
  numeric y() const;
  numeric z() const;

  numeric get_length() const;
  numeric get_length_squared() const;
  Vector unit() const;
  Vector vector_inverse() const;
  bool is_zero() const;
  Vector get_any_perpendicular() const;

  friend bool operator==(const Vector &a, const Vector &b) {
    auto diff = (a - b).get_length_squared();
    return diff < kEps * kEps;
  }
  friend bool operator!=(const Vector &a, const Vector &b) { return !(a == b); }
  friend Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a._x + b._x, a._y + b._y, a._z + b._z);
  }
  friend Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a._x - b._x, a._y - b._y, a._z - b._z);
  }
  friend numeric operator*(const Vector &a, const Vector &b)
  // dot product
  {
    return (a._x * b._x + a._y * b._y + a._z * b._z);
  }
  friend Vector operator^(const Vector &a, const Vector &b)
  // cross product
  {
    return Vector(a._y * b._z - a._z * b._y, a._z * b._x - a._x * b._z,
                  a._x * b._y - a._y * b._x);
  }
  friend Vector operator*(const Vector &a, const numeric k) {
    return Vector(a._x * k, a._y * k, a._z * k);
  }
  friend Vector operator*(const numeric k, const Vector &a) {
    return Vector(a._x * k, a._y * k, a._z * k);
  }
  friend Vector operator/(const Vector &a, const numeric k) {
    assertm(abs(k) >= kEps * kEps, "Division by 0.");
    numeric frac = inverse(k);
    return Vector(a._x * frac, a._y * frac, a._z * frac);
  }
  friend std::ostream &operator<<(std::ostream &os, const Vector &a) {
    os << '(' << a._x << ',' << a._y << ',' << a._z << ')';
    return os;
  }
};
