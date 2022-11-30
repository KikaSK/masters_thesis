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

  void _bound_consecutive(HalfEdge *previous, const HalfEdgeIndex i_previous,
                          HalfEdge *next, const HalfEdgeIndex i_next) const;
  void _bound_opposite(HalfEdge *AB, const HalfEdgeIndex i_AB, HalfEdge *BA,
                       const HalfEdgeIndex i_BA);
  void _bound_face(const FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                   HalfEdge *edge3) const;
  void _bound_opposite_outgoing(const MeshPoint &A, const MeshPointIndex i_B,
                                const HalfEdgeIndex i_BA);

 public:
  // TODO set as private after tests
  unordered_set<HalfEdge> _active_edges;
  unordered_set<HalfEdge> _checked_edges;
  unordered_set<HalfEdge> _mesh_edges_set;
  explicit Mesh(const Triangle &T);
  Mesh() = delete;
  Mesh(const Mesh &M) = delete;

  void cout_triangles() const;
  void add_triangle(HalfEdgeIndex index_AB, Point P, std::string type,
                    MeshPointIndex index_P = kInvalidPointIndex);

  bool is_active(HalfEdgeIndex index) const;
  bool is_checked(HalfEdgeIndex index) const;
  bool is_boundary(HalfEdgeIndex index) const;
  bool is_in_mesh(const HalfEdge &halfedge) const;

  HalfEdge get_halfedge(HalfEdgeIndex index) const;
  MeshPoint get_meshpoint(MeshPointIndex index) const;
  Face get_face(FaceIndex index) const;

  HalfEdge get_previous_halfedge(const HalfEdge &halfedge) const;
  HalfEdge get_previous_halfedge(HalfEdgeIndex index) const;
  HalfEdge get_next_halfedge(const HalfEdge &halfedge) const;
  HalfEdge get_next_halfedge(HalfEdgeIndex index) const;
  std::optional<HalfEdge> get_opposite_halfedge(const HalfEdge &halfedge) const;
  std::optional<HalfEdge> get_opposite_halfedge(HalfEdgeIndex index) const;

  HalfEdgeIndex get_previous_index(HalfEdgeIndex index) const;
  HalfEdgeIndex get_next_index(HalfEdgeIndex index) const;
  HalfEdgeIndex get_opposite_index(HalfEdgeIndex index) const;

  unordered_set<HalfEdge> get_active_edges() const;
  bool has_active_edge(const HalfEdge &halfedge) const;
  bool has_active_edge(const HalfEdgeIndex &index) const;
  size_t get_active_edges_size() const;

  void _assert_halfedge_index(const HalfEdgeIndex &index,
                              std::string function) const;
  void _assert_meshpoint_index(const MeshPointIndex &index,
                               std::string function) const;
  void _assert_face_index(const FaceIndex &index, std::string function) const;
};