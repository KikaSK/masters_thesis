#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H

#include <vector>

#include "algorithms.h"
#include "assertm.h"
#include "bounding_box.h"
#include "function.h"
#include "mesh.h"
#include "triangle.h"

class BasicAlgorithm {
 public:
  BasicAlgorithm(string name, Function f, Triangle seed_triangle,
                 numeric e_size, realsymbol x, realsymbol y, realsymbol z,
                 BoundingBox bounding_box)
      : name(name),
        F(f),
        my_mesh(seed_triangle),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
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

  Mesh calculate();
  void starting();
  void ending();
  bool step(const HalfEdgeIndex &working_edge);
  bool fix_proj(const HalfEdge &working_edge, const Point &projected,
                const Face &incident_face);
  /*bool fix_prev_next(const HalfEdge &working_edge, const Face &incident_face,
                     const MeshPoint &point, const bool is_prev,
                     const bool Delaunay);*/
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
  bool Delaunay_conditions(const HalfEdge &working_edge, const Point &P,
                           const Face &incident_face) const;
  bool non_Delaunay_conditions(const HalfEdge &working_edge, const Point &P,
                               const Face &incident_face) const;

  void create_triangle(const HalfEdge &working_edge, const Point &P,
                       const std::string type,
                       const MeshPointIndex index_P = -1);
  bool good_edges(const HalfEdge &working_edge, const Point &P) const;
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

  std::optional<Point> get_closest_point(const HalfEdge &working_edge,
                                         const Face &incident_face) const;
  std::optional<pair<Edge, numeric>> get_closest_edge(const Point &P,
                                                      const Triangle &N) const;

  void fix_corners();

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
  std::optional<HalfEdgeIndex> _find_closest_prev(
      const HalfEdge &working_edge, const Face &incident_face,
      const vector<HalfEdgeIndex> &prev) const;
  std::optional<HalfEdgeIndex> _find_closest_next(
      const HalfEdge &working_edge, const Face &incident_face,
      const vector<HalfEdgeIndex> &next) const;
};

#endif