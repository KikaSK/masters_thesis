#pragma once

#include "half_edge.h"
#include "triangle.h"

using HalfEdgeIndex = int;

class Face {
 private:
  Triangle _triangle;
  HalfEdgeIndex _halfedge;
  int _normal_multiplier;

 public:
  explicit Face(Triangle triangle, HalfEdgeIndex halfedge = kInvalidEdgeIndex);
  Face(const Face &F) = default;
  Face() = delete;

  Triangle get_triangle() const;
  HalfEdgeIndex get_halfedge() const;
  void set_halfedge(HalfEdgeIndex index);

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
};