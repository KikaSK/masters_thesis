#pragma once

#include "half_edge.h"
#include "triangle.h"

using HalfEdgeIndex = int;

class Face {
 private:
  HalfEdgeIndex _halfedge;
  Triangle _triangle;

 public:
  Face(HalfEdgeIndex halfedge, Triangle triangle);

  Triangle triangle() const;
  HalfEdgeIndex halfedge_index() const;
};