#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "face.h"
#include "half_edge.h"
#include "mesh_point.h"

using std::unordered_set;
using std::vector;

class Mesh {
 private:
  vector<MeshPoint> _mesh_points;
  vector<Face> _mesh_triangles;
  vector<HalfEdge> _mesh_edges;
  unordered_set<HalfEdge> _active_edges;
  unordered_set<HalfEdge> _checked_edges;
  unordered_set<HalfEdge> _mesh_edges_set;

  void _bound_consecutive(HalfEdge *previous, HalfEdgeIndex i_previous,
                          HalfEdge *next, HalfEdgeIndex i_next) const;
  void _bound_opposite(HalfEdge *AB, HalfEdgeIndex i_AB, HalfEdge *BA,
                       HalfEdgeIndex i_BA) const;
  void _bound_face(FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                   HalfEdge *edge3) const;
  void _bound_opposite_outgoing(const MeshPoint &A, MeshPointIndex i_B,
                                HalfEdge *BA, HalfEdgeIndex i_BA);

 public:
  Mesh(Face F);
  Mesh() = delete;

  void cout_triangles() const;
  void add_triangle(HalfEdgeIndex index_AB, Point P, std::string type,
                    MeshPointIndex index_P = kInvalidPointIndex);

  bool is_active(HalfEdgeIndex index) const;
  bool is_checked(HalfEdgeIndex index) const;
  bool is_boundary(HalfEdgeIndex index) const;
  bool is_active(HalfEdge halfedge) const;
  bool is_checked(HalfEdge halfedge) const;
  bool is_boundary(HalfEdge halfedge) const;
  bool is_in_mesh(HalfEdge halfedge) const;
};