#pragma once

#include <regex>
#include <vector>

#include "algorithms.h"
#include "assertm.h"
#include "bounding_box.h"
#include "function.h"
#include "mesh.h"
#include "singularity.h"
#include "triangle.h"

using GiNaC::realsymbol;
using std::pair;
using std::string;

class BasicAlgorithm {
 public:
  // the local mesh is created in the constructors
  explicit BasicAlgorithm(string name, Function f, Point seed_point,
                          numeric e_size, realsymbol x, realsymbol y,
                          realsymbol z, BoundingBox bounding_box,
                          const vector<Singularity> &singularities)
      : name(name),
        F(f),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
    std::cout << "Creating local mesh..." << endl;
    std::cout << "Number of singular points: " << singularities.size() << endl;
    for (int j = 0; j < singularities.size(); ++j) {
      int points_count = my_mesh.get_points_count();
      const Point singular_point = singularities[j].location();
      for (int i = 0; i < singularities[j].get_directions_count(); ++i) {
        // local mesh of singular point
        triangulate_singular_point_local(
            singularities[j], i, (i == 0) ? kInvalidPointIndex : points_count);
      }
    }
  }

  explicit BasicAlgorithm(string name, Function f, const Function &F,
                          const Function &G,
                          const vector<vector<Point>> &polylines,
                          numeric e_size, realsymbol x, realsymbol y,
                          realsymbol z, BoundingBox bounding_box,
                          const string type)
      : name(name),
        F(f),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
    for (auto polyline : polylines)
        // local mesh of singular curve
      get_local_mesh_curve(F, G, f, polyline, bounding_box, e_size, type);
  }

  explicit BasicAlgorithm(string name, Function f, Triangle seed_triangle,
                          numeric e_size, realsymbol x, realsymbol y,
                          realsymbol z, BoundingBox bounding_box)
      : name(name),
        F(f),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
    // first triangle - regular surface starting
    my_mesh.add_first_triangle(seed_triangle, bounding_box);
  }

  void calculate();
  void starting();
  void ending();
  bool step(const HalfEdgeIndex &working_edge);
  bool fix_proj(const HalfEdge &working_edge, const Point &projected,
                const bool Delaunay);
  bool fix_prev_next(const HalfEdge &working_edge, const bool Delaunay);
  bool fix_close_points(const HalfEdge &working_edge, const Point &projected,
                        const bool Delaunay);
  bool fix_overlap(const HalfEdge &working_edge, const MeshPoint &overlap_point,
                   const bool Delaunay);
  bool fix_holes(const HalfEdge &working_edge);
  bool fix_breakers(const HalfEdge &working_edge, const Point &projected,
                    const bool Delaunay);

  MeshPointIndex find_triangle_hole(const HalfEdge &working_edge);
  std::pair<FaceIndex, FaceIndex> find_neighbor_faces(
      const HalfEdge &working_edge, const Point &P) const;
  std::optional<vector<MeshPoint>> find_close_points(
      const Point &P, const HalfEdge &working_edge) const;
  bool basic_triangle(const HalfEdge &working_edge, const MeshPoint &point);

  // conditions
  bool check_conditions(const HalfEdge &working_edge, const Point &P,
                        const bool Delaunay,
                        const MeshPointIndex i_P = kInvalidPointIndex) const;
  bool Delaunay_conditions(const HalfEdge &working_edge, const Point &P) const;
  bool non_Delaunay_conditions(const HalfEdge &working_edge, const Point &P,
                               const MeshPointIndex i_P) const;

  bool neighbor_triangles_normal_check(const Triangle &T1,
                                       const Triangle &T2) const;
  bool orientability_check(const HalfEdge &working_edge,
                           const MeshPoint &point) const;
  bool overlap_normals_check(const MeshPoint &candidate,
                             const HalfEdge &working_edge) const;
  bool gc_distance_check(const Triangle T) const;
  bool good_edges(const HalfEdge &working_edge, const Point &P) const;
  bool good_new_point(const MeshPoint &point) const;

  std::optional<Point> get_projected(const HalfEdge &working_edge) const;
  void create_triangle(const HalfEdge &working_edge, const Point &P,
                       const bool is_new, const MeshPointIndex index_P = -1);

  bool is_active(const HalfEdgeIndex &halfedge_index) const;
  bool is_checked(const HalfEdgeIndex &halfedge_index) const;
  bool is_bounding(const HalfEdgeIndex &halfedge_index) const;
  bool is_border(const HalfEdgeIndex &halfedge_index) const;
  void delete_from_active(const HalfEdgeIndex &halfedge_index);
  void delete_from_checked(const HalfEdgeIndex &halfedge_index);
  void push_edge_to_active(const HalfEdgeIndex &halfedge_index);
  void push_edge_to_checked(const HalfEdgeIndex &halfedge_index);

  //void fix_corners();

  // singularities
  void triangulate_singular_point_local(const Singularity &singularity,
                                        const int branch,
                                        const MeshPointIndex singular_index);
  void get_points_singular_point(const Singularity singularity,
                                 const int num_triangles,
                                 const numeric &edge_length,
                                 vector<vector<Point>> &points,
                                 const int branch = 0, const int layers = 1);
  void get_local_mesh_point(const vector<vector<Point>> &points,
                            const Singularity &singularity,
                            const MeshPointIndex singular_index,
                            const int num_triangles, const int branch,
                            const int layers);
  void get_local_mesh_curve(const Function &_F, const Function &_G,
                            const Function &inter_FG,
                            const vector<Point> &polyline,
                            const BoundingBox &bounding_box,
                            const numeric e_size, const string type);
  
 private:
  string name;
  const Function F;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
  BoundingBox bounding_box;
  bool _fix_prev_next(const HalfEdge &working_edge, const MeshPoint &prev,
                      const bool Delaunay);
  std::optional<vector<MeshPoint>> _linear_close_points_finder(
      const Point &P, const HalfEdge &working_edge) const;
  std::optional<vector<MeshPoint>> _tree_close_points_finder(
      const Point &P, const HalfEdge &working_edge) const;
  bool _check_normals(const HalfEdge &working_edge, const Point &point) const;
};
