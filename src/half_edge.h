#pragma once

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
  HalfEdgeIndex _self;
  bool _is_active;
  bool _is_checked;

 public:
  HalfEdge(Edge edge, MeshPointIndex A, MeshPointIndex B,
           HalfEdgeIndex previous = kInvalidEdgeIndex,
           HalfEdgeIndex next = kInvalidEdgeIndex,
           HalfEdgeIndex opposite = kInvalidEdgeIndex,
           FaceIndex incident = kInvalidFaceIndex);
  HalfEdge() = delete;

  MeshPointIndex get_A() const;
  MeshPointIndex get_B() const;

  HalfEdgeIndex get_opposite() const;
  void set_opposite(HalfEdgeIndex opposite);
  HalfEdgeIndex get_previous() const;
  void set_previous(HalfEdgeIndex previous);
  HalfEdgeIndex get_next() const;
  void set_next(HalfEdgeIndex next);
  FaceIndex get_incident() const;
  void set_incident(FaceIndex incident);
  HalfEdgeIndex get_index() const;
  void set_index(HalfEdgeIndex index);

  void set_active();
  void set_checked();
  void set_inside();

  bool is_active() const;
  bool is_checked() const;
  bool is_boundary() const;
  bool is_inside() const;

  numeric get_length() const;
  Point get_midpoint() const;

  Edge get_edge() const;
  Point get_point_A() const;
  Point get_point_B() const;

  // partial equality of half edges, opposite edges are partially equal
  friend bool operator%(const HalfEdge &e1, const HalfEdge &e2) {
    return (e1._edge == e2._edge);
  }
  // total equality of half edges, including the direction
  friend bool operator==(const HalfEdge &e1, const HalfEdge &e2) {
    return (e1._A == e2._A && e1._B == e2._B);
  }
  friend bool operator!=(const HalfEdge &e1, const HalfEdge &e2) {
    return !(e1 == e2);
  }

  friend std::ostream &operator<<(std::ostream &os, const HalfEdge &e) {
    os << "HalfEdge AB: " << endl;
    os << "A: " << e.get_point_A() << endl << "B: " << e.get_point_B() << endl;
    os << "Previous: " << e.get_previous() << " Next: " << e.get_next()
       << " Opposite: " << e.get_opposite() << endl;
    return os;
  }
};

namespace std {
template <>
struct hash<HalfEdge> {
  std::size_t operator()(const HalfEdge &e) const {
    using std::hash;
    /*
    HalfEdgeIndex opposite = e.get_opposite();
    HalfEdgeIndex previous = e.get_previous();
    HalfEdgeIndex next = e.get_next();
    FaceIndex incident = e.get_incident();
    return hash<HalfEdgeIndex>()(opposite) ^ hash<HalfEdgeIndex>()(previous) ^
           hash<HalfEdgeIndex>()(next) ^ hash<HalfEdgeIndex>()(incident);
    */
    int Ax = (e.get_edge().A().x() * pow(10, 6)).to_int();
    int Ay = (e.get_edge().A().y() * pow(10, 6)).to_int();
    int Az = (e.get_edge().A().z() * pow(10, 6)).to_int();
    int Bx = (e.get_edge().B().x() * pow(10, 6)).to_int();
    int By = (e.get_edge().B().y() * pow(10, 6)).to_int();
    int Bz = (e.get_edge().B().z() * pow(10, 6)).to_int();
    return hash<int>()(Ax) ^ hash<int>()(Ay) ^ hash<int>()(Az) ^
           hash<int>()(Bx) ^ hash<int>()(By) ^ hash<int>()(Bz);
  }
};
}  // namespace std
