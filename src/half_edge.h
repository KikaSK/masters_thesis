#ifndef HALF_EDGE_H
#define HALF_EDGE_H

#include <ginac/ginac.h>

#include <iostream>
#include <optional>

#include "assertm.h"
#include "constants.h"
#include "edge.h"
#include "mesh_point.h"
#include "point.h"
#include "vector.h"

using namespace GiNaC;
using std::endl;

using MeshPointIndex = int;
using HalfEdgeIndex = int;
using MeshTriangleIndex = int;
using FaceIndex = int;

// class Triangle;
class Point;
class MeshPoint;
class Vector;
class Edge;

class HalfEdge {
 private:
  Edge _edge;
  MeshPointIndex _A;
  // todo remove point B
  MeshPointIndex _B;
  HalfEdgeIndex _previous;
  HalfEdgeIndex _next;
  HalfEdgeIndex _opposite;
  FaceIndex _incident;

 public:
  HalfEdge(Edge edge, MeshPointIndex A, MeshPointIndex B,
           HalfEdgeIndex previous = -1, HalfEdgeIndex next = -1,
           HalfEdgeIndex opposite = -1, FaceIndex incident = -1);
  HalfEdge() = delete;

  MeshPointIndex A() const;
  MeshPointIndex B() const;

  HalfEdgeIndex get_opposite() const;
  void set_opposite(HalfEdgeIndex opposite);
  HalfEdgeIndex get_previous() const;
  void set_previous(HalfEdgeIndex previous);
  HalfEdgeIndex get_next() const;
  void set_next(HalfEdgeIndex next);
  FaceIndex get_incident() const;
  void set_incident(FaceIndex incident);

  numeric get_length() const;
  Point get_midpoint() const;

  friend bool operator==(const HalfEdge &e1, const HalfEdge &e2) {
    return (e1._A == e2._A && e1._B == e2._B) ||
           (e1._A == e2._B && e1._B == e2._A);
  }
  friend bool operator!=(const HalfEdge &e1, const HalfEdge &e2) {
    return !(e1 == e2);
  }

  friend std::ostream &operator<<(std::ostream &os, const HalfEdge &e) {
    os << "A: " << e._A << endl << "B: " << e._B << endl;
    return os;
  }
};

#endif
