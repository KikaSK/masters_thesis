#ifndef MESH_H
#define MESH_H

#include <vector>

#include "face.h"
#include "half_edge.h"
#include "mesh_point.h"

using std::vector;

class Mesh {
 private:
  vector<MeshPoint> _mesh_points;
  vector<Face> _mesh_triangles;
  vector<HalfEdge> _mesh_edges;

 public:
  Mesh(Face F);
  Mesh() = delete;

  void cout_triangles() const;
  void add_triangle(Edge e, Point P);
};

#endif