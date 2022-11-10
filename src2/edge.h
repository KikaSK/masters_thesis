#ifndef EDGE_H
#define EDGE_H

#include <ginac/ginac.h>

#include <iostream>
#include <optional>

#include "assertm.h"
#include "mesh_point.h"
#include "point.h"
#include "vector.h"

using namespace GiNaC;
using std::endl;

// class Triangle;
class Point;
class MeshPoint;
class Vector;

class Edge {
 private:
  MeshPoint *_A;
  // todo remove point B
  MeshPoint *_B;
  Edge *_previous;
  Edge *_next;
  Edge *_opposite;
  // Face* _incident;

 public:
  Edge(MeshPoint *A, Vector u, Edge *previous = nullptr, Edge *next = nullptr,
       Edge *opposite = nullptr);
  Edge(MeshPoint *A, MeshPoint *B, Edge *previous = nullptr,
       Edge *next = nullptr, Edge *opposite = nullptr);
  Edge() = delete;

  MeshPoint *A() const;
  MeshPoint *B() const;

  Edge *get_opposite() const;
  void set_opposite(Edge *const opposite);
  Edge *get_previous() const;
  void set_previous(Edge *const previous);
  Edge *get_next() const;
  void set_next(Edge *const next);

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
