#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <ginac/ginac.h>

#include <cmath>
#include <iostream>

#include "assertm.h"
#include "constants.h"
#include "edge.h"
#include "point.h"

using namespace GiNaC;
using std::endl;

class Triangle {
 private:
  Point _A, _B, _C;
  Edge _AB, _BC, _CA;

 public:
  Triangle(Point A, Point B, Point C);
  Triangle() = delete;

  Point A() const;
  Point B() const;
  Point C() const;

  Edge AB() const;
  Edge BC() const;
  Edge CA() const;

  Point get_gravity_center() const;
  bool is_triangle() const;
  Point get_circumcenter() const;
  Vector get_normal() const;
  bool is_in_triangle(Point P) const;

  // a b c
  // a c b
  // b a c
  // b c a
  // c a b
  // c b a

  friend bool operator==(const Triangle &a, const Triangle &b) {
    return ((a._A == b._A && a._B == b._B && a._C == b._C) ||
            (a._A == b._A && a._B == b._C && a._C == b._B) ||
            (a._A == b._B && a._B == b._A && a._C == b._C) ||
            (a._A == b._B && a._B == b._C && a._C == b._A) ||
            (a._A == b._C && a._B == b._A && a._C == b._B) ||
            (a._A == b._C && a._B == b._B && a._C == b._A));
  }
  friend std::ostream &operator<<(std::ostream &os, const Triangle &T) {
    os << "A: " << T._A << endl
       << "B: " << T._B << endl
       << "C: " << T._C << endl;
    return os;
  }
};

#endif