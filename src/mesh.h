#pragma once

#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "../range-tree/nd_reb.h"
#include "bounding_box.h"
#include "face.h"
#include "half_edge.h"
#include "mesh_point.h"

using std::unordered_set;
using std::vector;

enum class NewTriangleType {
  kNew = 0,
  kPrevious = 1,
  kNext = 2,
  kOverlap = 3,
  kFill = 4,
};

class Mesh {
 private:
  vector<MeshPoint> _mesh_points;
  vector<Face> _mesh_triangles;
  vector<HalfEdge> _mesh_edges;
  unordered_set<HalfEdgeIndex> _active_edges;
  unordered_set<HalfEdgeIndex> _checked_edges;
  unordered_set<HalfEdgeIndex> _mesh_edges_set;
  unordered_set<HalfEdgeIndex> _bounding_edges;
  Tree<MeshPoint> _point_tree;

  void _bound_consecutive(HalfEdge *previous, const HalfEdgeIndex i_previous,
                          HalfEdge *next, const HalfEdgeIndex i_next) const;
  void _bound_opposite(HalfEdge *AB, const HalfEdgeIndex i_AB, HalfEdge *BA,
                       const HalfEdgeIndex i_BA);
  void _bound_face(const FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                   HalfEdge *edge3) const;
  void _bound_opposite_outgoing(const MeshPoint &A, const MeshPointIndex i_B,
                                const HalfEdgeIndex i_BA);
  bool _tree_is_in_mesh(const Edge &edge) const;
  bool _linear_is_in_mesh(const Edge &edge) const;
  HalfEdgeIndex _tree_get_edge_index(const Edge &edge) const;
  HalfEdgeIndex _linear_get_edge_index(const Edge &edge) const;

 public:
  explicit Mesh(const Triangle &T, const BoundingBox &bounding_box);
  Mesh() = default;
  Mesh(const Mesh &M) = delete;
  Mesh(Mesh &&M) = default;
  Mesh &operator=(Mesh &&M) = default;

  void cout_triangles() const;
  void cout_mesh() const;
  void cout_triangles_number() const;

  void add_first_triangle(const Triangle &T, const BoundingBox &bounding_box);
  void add_triangle_to_meshpoint(MeshPointIndex i_A, Point point_B,
                                 Point point_C,
                                 const BoundingBox &bounding_box);
  void add_triangle_to_two_meshpoints(MeshPointIndex i_A, MeshPointIndex i_B,
                                      Point point_C,
                                      const BoundingBox &bounding_box);
  void add_triangle(HalfEdgeIndex index_AB, Point P, const bool is_new,
                    const BoundingBox &bounding_box,
                    MeshPointIndex index_P = kInvalidPointIndex);

  bool is_active(HalfEdgeIndex index) const;
  bool is_checked(HalfEdgeIndex index) const;
  bool is_bounding(HalfEdgeIndex index) const;
  bool is_boundary(HalfEdgeIndex index) const;
  bool is_in_mesh(const Edge &edge) const;
  bool is_boundary_point(const MeshPoint &P) const;

  HalfEdge get_halfedge(HalfEdgeIndex index) const;
  MeshPoint get_meshpoint(MeshPointIndex index) const;
  Face get_face(FaceIndex index) const;

  int get_faces_count() const;
  int get_edges_count() const;
  int get_points_count() const;

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
  HalfEdgeIndex get_edge_index(const Edge &edge) const;
  vector<MeshPoint> get_prev(const HalfEdge &halfedge) const;
  vector<MeshPoint> get_next(const HalfEdge &halfedge) const;

  unordered_set<HalfEdgeIndex> get_active_edges() const;
  void set_active_edges(const unordered_set<HalfEdgeIndex> &edges);
  unordered_set<HalfEdgeIndex> get_checked_edges() const;
  void set_checked_edges(const unordered_set<HalfEdgeIndex> &edges);
  unordered_set<HalfEdgeIndex> get_bounding_edges() const;
  void move_active_edges_to_checked();
  vector<MeshPoint> get_mesh_points() const;
  vector<HalfEdge> get_mesh_edges() const;
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
  bool has_outgoing_next(const HalfEdge &working_edge,
                         const Point &point) const;
  bool has_incoming_prev(const HalfEdge &working_edge,
                         const Point &point) const;

  size_t get_active_edges_size() const;
  size_t get_checked_edges_size() const;
  bool active_edges_empty() const;
  bool checked_edges_empty() const;

  bool check_Delaunay(const HalfEdge &working_edge,
                      const Point &new_point) const;
  vector<MeshPoint> get_breakers(const Triangle &T) const;
  vector<MeshPoint> get_meshpoints_in_interval(numeric min_x, numeric max_x,
                                               numeric min_y, numeric max_y,
                                               numeric min_z,
                                               numeric max_z) const;

  void obj_format(const std::string &name) const;

  bool edges_check(const std::string &message,
                   const HalfEdgeIndex working_edge = kInvalidEdgeIndex) const;
  NewTriangleType _find_type(const HalfEdgeIndex index_AB,
                             const MeshPoint &P) const;
  vector<MeshPoint> _linear_breakers_getter(const Triangle &T) const;
  vector<MeshPoint> _tree_breakers_getter(const Triangle &T) const;

  // moved from algorithm.h

  // angle BAP in range (-Pi, Pi) with respect to neighbour triangle
  std::optional<numeric> angle(const HalfEdge &working_edge, const Point &P,
                               const bool clockwise) const;

  // true if angle is between 0 and 3*pi/4 with respect to neighbour triangle
  bool good_orientation(const HalfEdge &working_edge, const Point &P) const;

  numeric _midpoint_line_point_distance(const HalfEdge &working_edge,
                                        const Point &P) const;
  numeric _linesegment_line_point_distance(const HalfEdge &working_edge,
                                           const Point &P) const;
  // https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
  // returns ditance between point and line given by working edge
  numeric line_point_dist(const HalfEdge &working_edge, const Point &P) const;
};