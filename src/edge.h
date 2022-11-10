#ifndef EDGE_H
#define EDGE_H

#include <ginac/ginac.h>

#include <iostream>
#include <optional>

#include "assertm.h"
#include "constants.h"
#include "point.h"
#include "vector.h"

using namespace GiNaC;
using std::endl;

// class Triangle;
class Point;
class Vector;

class Edge {
 private:
  Point _A;
  Point _B;

 public:
  Edge(Point A, Vector u);
  Edge(Point A, Point B);
  Edge() = delete;

  Point A() const;
  Point B() const;

  numeric get_length() const;
  Point get_midpoint() const;

  friend bool operator==(const Edge &e1, const Edge &e2) {
    return (e1._A == e2._A && e1._B == e2._B) ||
           (e1._A == e2._B && e1._B == e2._A);
  }
  friend bool operator!=(const Edge &e1, const Edge &e2) { return !(e1 == e2); }

  friend std::ostream &operator<<(std::ostream &os, const Edge &e) {
    os << "A: " << e._A << endl << "B: " << e._B << endl;
    return os;
  }
};

#endif
