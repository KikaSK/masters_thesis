#pragma once

#include <vector>

#include "algorithms.h"
#include "assertm.h"
#include "bounding_box.h"
#include "function.h"
#include "mesh.h"
#include "triangle.h"

using GiNaC::realsymbol;
using std::pair;
using std::string;

class BasicAlgorithm {
 public:
  BasicAlgorithm(string name, Function f, Point seed_point, numeric e_size,
                 realsymbol x, realsymbol y, realsymbol z,
                 BoundingBox bounding_box)
      : name(name),
        F(f),
        // my_mesh(seed_triangle),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
    Triangle seed_triangle = find_seed_triangle(seed_point);
    std::cout << "seed triangle created" << endl;
    my_mesh = Mesh(seed_triangle, bounding_box);
    std::cout << "in basic algorithm constructor" << endl;
    /*
    if (bounding_box.new_bounding_edge(seed_triangle.AB()))
      bounding_edges.push_back(seed_triangle.AB());
    else
      active_edges.push_back(seed_triangle.AB());

    if (bounding_box.new_bounding_edge(seed_triangle.BC()))
      bounding_edges.push_back(seed_triangle.BC());
    else
      active_edges.push_back(seed_triangle.BC());

    if (bounding_box.new_bounding_edge(seed_triangle.CA()))
      bounding_edges.push_back(seed_triangle.CA());
    else
      active_edges.push_back(seed_triangle.CA());
    */
  }

  void calculate();
  void starting();
  void ending();
  bool step(const HalfEdgeIndex &working_edge);
  bool fix_proj(const HalfEdge &working_edge, const Point &projected,
                const Face &incident_face, const bool Delaunay);
  bool fix_prev_next(const HalfEdge &working_edge, const bool Delaunay);
  bool fix_close_points(const HalfEdge &working_edge, const Face &incident_face,
                        const Point &projected);
  bool fix_overlap(const HalfEdge &working_edge, const Face &incident_face,
                   const MeshPoint &overlap_point, const bool Delaunay);
  int fix_holes(const HalfEdge &working_edge, const Face &incident_face);
  bool fix_holes2(const HalfEdge &working_edge);
  bool fix_breakers(const HalfEdge &working_edge, const Point &projected,
                    const Face &incident_face, const bool Delaunay);
  bool fix_same_points(const HalfEdge &working_edge, const Point &projected,
                       const Face &incident_face);

  pair<MeshPoint, MeshPoint> find_prev_next(const HalfEdge &working_edge,
                                            const Face &incident_face) const;
  MeshPointIndex find_triangle_hole(const HalfEdge &working_edge);
  std::pair<FaceIndex, FaceIndex> find_neighbor_faces(
      const HalfEdge &working_edge, const MeshPoint &P) const;
  bool neighbor_triangles_normal_check(const Triangle &T1,
                                       const Triangle &T2) const;
  bool orientability_check(const HalfEdge &working_edge,
                           const MeshPoint &point) const;

  Point get_projected(const HalfEdge &working_edge,
                      const Face &incident_face) const;
  bool overlap_normals_check(const MeshPoint &candidate,
                             const HalfEdge &working_edge) const;
  std::optional<vector<MeshPoint>> find_close_points(
      const Point &P, const HalfEdge &working_edge,
      const Face &incident_face) const;
  int update_border(const Edge &new_edge1, const Edge &new_edge2);
  bool check_conditions(const HalfEdge &working_edge, const Point &P,
                        const bool Delaunay) const;
  bool Delaunay_conditions(const HalfEdge &working_edge, const Point &P) const;
  bool non_Delaunay_conditions(const HalfEdge &working_edge,
                               const Point &P) const;

  void create_triangle(const HalfEdge &working_edge, const Point &P,
                       const bool is_new, const MeshPointIndex index_P = -1);
  bool good_edges(const HalfEdge &working_edge, const Point &P) const;
  bool good_new_point(const MeshPoint &point) const;
  bool basic_triangle(const HalfEdge &working_edge, const Face &incident_face,
                      const MeshPoint &point);

  bool is_active(const HalfEdgeIndex &halfedge_index) const;
  bool is_checked(const HalfEdgeIndex &halfedge_index) const;
  bool is_bounding(const HalfEdgeIndex &halfedge_index) const;
  bool is_border(const HalfEdgeIndex &halfedge_index) const;
  bool is_border_point(const MeshPoint &P) const;
  bool is_border_point(const Point &P) const;
  void delete_from_active(const HalfEdgeIndex &halfedge_index);
  void delete_from_checked(const HalfEdgeIndex &halfedge_index);
  void push_edge_to_active(const HalfEdgeIndex &halfedge_index);
  void push_edge_to_checked(const HalfEdgeIndex &halfedge_index);

  std::optional<MeshPoint> get_closest_point(const HalfEdge &working_edge,
                                             const Face &incident_face) const;
  std::optional<pair<Edge, numeric>> get_closest_edge(const Point &P,
                                                      const Triangle &N) const;

  void fix_corners();

  // section seed triangle
  Edge get_seed_edge(Point seed_point) const;
  Point get_seed_triangle(const Edge &e) const;
  Triangle find_seed_triangle(Point seed) const;
  void create_mesh(Point seed);

 private:
  string name;
  Function F;
  // vector<Edge> active_edges;
  // vector<Edge> checked_edges;
  // vector<Edge> bounding_edges;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
  BoundingBox bounding_box;
  HalfEdgeIndex _find_prev(const HalfEdge &working_edge,
                           const Face &incident_face) const;
  HalfEdgeIndex _find_next(const HalfEdge &working_edge,
                           const Face &incident_face) const;
  bool _fix_prev_next(const HalfEdge &working_edge, const MeshPoint &prev,
                      const bool Delaunay);
  std::optional<HalfEdgeIndex> _find_closest_prev(
      const HalfEdge &working_edge, const Face &incident_face,
      const vector<HalfEdgeIndex> &prev) const;
  std::optional<HalfEdgeIndex> _find_closest_next(
      const HalfEdge &working_edge, const Face &incident_face,
      const vector<HalfEdgeIndex> &next) const;
  std::optional<vector<MeshPoint>> _linear_close_points_finder(
      const Point &P, const HalfEdge &working_edge,
      const Face &incident_face) const;
  std::optional<vector<MeshPoint>> _tree_close_points_finder(
      const Point &P, const HalfEdge &working_edge,
      const Face &incident_face) const;
  bool _check_normals(const HalfEdge &working_edge, const Point &point) const;
};
