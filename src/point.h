#pragma once

#include <ginac/ginac.h>

#include <cmath>
#include <iostream>

#include "assertm.h"
#include "constants.h"
#include "vector.h"

using GiNaC::numeric;
using std::endl;

class Vector;

class Point {
 private:
  numeric _x, _y, _z;

  // TODO Ajko rule of 3

 public:
  Point(const numeric x, const numeric y, const numeric z);
  Point(const Point &A, const Vector &u);
  // default constructor
  Point() = default;
  Point(const Point &P) = default;

  numeric x() const;
  numeric y() const;
  numeric z() const;

  friend std::ostream &operator<<(std::ostream &os, const Point &A) {
    os << '[' << (A._x).to_double() << ',' << (A._y).to_double() << ','
       << (A._z).to_double() << ']';
    return os;
  }
  friend std::ostream &operator<<(std::ostream &os, const Point *const A) {
    os << '[' << A->_x << ',' << A->_y << ',' << A->_z << ']';
    return os;
  }
  friend Point operator+(const Point &A, const Point &B) {
    return Point(A._x + B._x, A._y + B._y, A._z + B._z);
  }
  friend bool operator==(const Point &A, const Point &B) {
    auto diff = GiNaC::pow(B._x - A._x, 2) + GiNaC::pow(B._y - A._y, 2) +
                GiNaC::pow(B._z - A._z, 2);
    return diff < kEps * kEps;
  }
  friend bool operator!=(const Point &A, const Point &B) { return !(A == B); }
  friend Point operator*(const numeric k, const Point &A) {
    return Point(A._x * k, A._y * k, A._z * k);
  }
};
