#pragma once

#include "half_edge.h"
#include "triangle.h"

using HalfEdgeIndex = int;

class Face {
 private:
  Triangle _triangle;
  HalfEdgeIndex _halfedge;

 public:
  Face(Triangle triangle, HalfEdgeIndex halfedge = kInvalidEdgeIndex);
  Face(const Face &F) = default;
  Face() = delete;

  Triangle get_triangle() const;
  HalfEdgeIndex get_halfedge() const;
  void set_halfedge(HalfEdgeIndex index);
};