#pragma once

#include <fstream>
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
  unordered_set<HalfEdgeIndex> _active_edges;
  unordered_set<HalfEdgeIndex> _checked_edges;
  unordered_set<HalfEdgeIndex> _mesh_edges_set;
  unordered_set<HalfEdgeIndex> _bounding_edges;

  void _bound_consecutive(HalfEdge *previous, const HalfEdgeIndex i_previous,
                          HalfEdge *next, const HalfEdgeIndex i_next) const;
  void _bound_opposite(HalfEdge *AB, const HalfEdgeIndex i_AB, HalfEdge *BA,
                       const HalfEdgeIndex i_BA);
  void _bound_face(const FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                   HalfEdge *edge3) const;
  void _bound_opposite_outgoing(const MeshPoint &A, const MeshPointIndex i_B,
                                const HalfEdgeIndex i_BA);

 public:
  explicit Mesh(const Triangle &T);
  Mesh() = delete;
  Mesh(const Mesh &M) = delete;
  Mesh(Mesh &&M) = default;

  void cout_triangles() const;
  void cout_mesh() const;
  void cout_triangles_number() const;

  void add_triangle(HalfEdgeIndex index_AB, Point P, std::string type,
                    MeshPointIndex index_P = kInvalidPointIndex);

  bool is_active(HalfEdgeIndex index) const;
  bool is_checked(HalfEdgeIndex index) const;
  bool is_bounding(HalfEdgeIndex index) const;
  bool is_boundary(HalfEdgeIndex index) const;
  bool is_in_mesh(const Edge &edge) const;
  HalfEdgeIndex get_edge_index(const Edge &edge) const;
  bool is_boundary_point(const MeshPoint &P) const;

  HalfEdge get_halfedge(HalfEdgeIndex index) const;
  MeshPoint get_meshpoint(MeshPointIndex index) const;
  Face get_face(FaceIndex index) const;

  HalfEdge get_previous_halfedge(const HalfEdge &halfedge) const;
  HalfEdge get_previous_halfedge(HalfEdgeIndex index) const;
  HalfEdge get_next_halfedge(const HalfEdge &halfedge) const;
  HalfEdge get_next_halfedge(HalfEdgeIndex index) const;
  HalfEdgeIndex get_previous_index(HalfEdgeIndex index) const;
  HalfEdgeIndex get_next_index(HalfEdgeIndex index) const;
  HalfEdgeIndex get_opposite_index(HalfEdgeIndex index) const;
  std::optional<HalfEdge> get_opposite_halfedge(const HalfEdge &halfedge) const;
  std::optional<HalfEdge> get_opposite_halfedge(HalfEdgeIndex index) const;
  FaceIndex get_incident_face(const HalfEdgeIndex &index) const;

  unordered_set<HalfEdgeIndex> get_active_edges() const;
  unordered_set<HalfEdgeIndex> get_checked_edges() const;
  unordered_set<HalfEdgeIndex> get_bounding_edges() const;
  vector<MeshPoint> get_mesh_points() const;
  vector<Face> get_mesh_faces() const;
  HalfEdgeIndex get_active_edge() const;
  HalfEdgeIndex get_checked_edge() const;
  void remove_active_edge(const HalfEdgeIndex &index);
  void add_edge_to_active(const HalfEdgeIndex &index);
  void remove_checked_edge(const HalfEdgeIndex &index);
  void add_edge_to_checked(const HalfEdgeIndex &index);

  bool has_active_edge(const HalfEdge &halfedge) const;
  bool has_active_edge(const HalfEdgeIndex &index) const;
  bool has_checked_edge(const HalfEdge &halfedge) const;
  bool has_checked_edge(const HalfEdgeIndex &index) const;

  size_t get_active_edges_size() const;
  size_t get_checked_edges_size() const;
  bool active_edges_empty() const;
  bool checked_edges_empty() const;

  bool check_Delaunay(const Mesh &mesh, const Triangle &new_triangle,
                      const HalfEdge &working_edge,
                      const Face &incident_face) const;

  void obj_format(const std::string &name) const;
};