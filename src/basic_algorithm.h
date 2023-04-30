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
    // Vector singular_direction(0, 0, 1);
    // Vector singular_direction(-1, 0.5, 0);
    // const Point singular_point(0, 0, 0);
    std::cout << "Creating local mesh: " << endl;
    std::cout << "Number of singular points: " << singularities.size() << endl;

    std::regex A1(".*A[0-9]*[02468][+][-].*");
    std::regex D1(".*D[0-9]*[02468][+][-].*");
    std::regex E1(".*((E6[+][+])|(E8[+][+])|(E6[+]-)).*");

    std::regex A2(".*A[0-9]*[13579][+][-].*");
    std::regex D2(".*D[0-9]*[13579][+][-].*");

    std::regex A3(".*A[0-9]+[-][-].*");

    for (int j = 0; j < singularities.size(); ++j) {
      const Point singular_point = singularities[j].location();
      for (int i = 0; i < singularities[j].get_directions_count(); ++i) {
        Vector singular_direction = singularities[j].get_direction(i);
        std::cout << "Sing dir: " << singular_direction << endl;
        if (std::regex_match(name, A1) || 
            std::regex_match(name, D1) ||
            std::regex_match(name, E1)  /*||
            name == "./outputs/my_run_input_A2+-_0.9" ||
            name == "./outputs/TEXT_A2+-_0.5" ||
            name == "./outputs/TEXT_A2-example-1_0.2" ||
            name == "./outputs/my_run_input_A2-example-1_0.2" ||
            name == "./outputs/my_run_input_A4+-_1"*/) {
          std::cout << "group1" << endl;

          triangulate_singularity_case2(singularities[j], i,
                                        (i == 0) ? kInvalidPointIndex : 0);
        } else if (std::regex_match(name, A2) /*||
std::regex_match(name, D2)
name == "./outputs/my_run_input_D4--_0.4" ||
name == "./outputs/my_run_input_D4--_0.2" ||
name == "./outputs/my_run_input_D4--_0.5" ||
name == "./outputs/my_run_input_D4--_0.6" ||
name == "./outputs/my_run_input_D4--_0.1" ||
name == "./outputs/my_run_input_D4--_0.3"
*/) {
          std::cout << "group2" << endl;
          triangulate_singularity_circular(singularities[j], i,
                                           (i == 0) ? kInvalidPointIndex : 0);
        } else if (std::regex_match(name, A3)) {
          std::cout << "group3" << endl;
          triangulate_An_analytical(singularities[j], i,
                                    (i == 0) ? kInvalidPointIndex : 0);
        } else {
          std::cout << "group4" << endl;
          triangulate_cone_iterative(singular_point, singular_direction,
                                     &my_mesh,
                                     (i == 0) ? kInvalidPointIndex : 0);
        }
      }
    }
  }
  /*
  explicit BasicAlgorithm(string name, Function f, Point seed_point,
                          numeric e_size, realsymbol x, realsymbol y,
                          realsymbol z, BoundingBox bounding_box,
                          const vector<Point> &singular_points,
                          const vector<vector<Vector>> &singular_directions,
                          const vector<int> &types = {})
      : name(name),
        F(f),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
    // Vector singular_direction(0, 0, 1);
    // Vector singular_direction(-1, 0.5, 0);
    // const Point singular_point(0, 0, 0);
    std::cout << "Types: " << types.empty() << endl;
    std::cout << "Creating local mesh: " << endl;
    std::cout << "Number of singular points: " << singular_directions.size()
              << endl;

    std::regex A1(".*A[0-9]*[02468][+][-].*");
    std::regex D1(".*D[0-9]*[02468][+][-].*");
    std::regex E1(".*((E6[+][+])|(E8[+][+])|(E6[+]-)).*");

    std::regex A2(".*A[0-9]*[13579][+][-].*");
    std::regex D2(".*D[0-9]*[13579][+][-].*");

    std::regex A3(".*A[0-9]+[-][-].*");

    for (int j = 0; j < singular_directions.size(); ++j) {
      const Point singular_point = singular_points[j];
      for (int i = 0; i < singular_directions[j].size(); ++i) {
        Vector singular_direction = singular_directions[j][i];
        std::cout << "Sing dir: " << singular_direction << endl;
        if (std::regex_match(name, A1) ||
            std::regex_match(name, D1) ||
            std::regex_match(name, E1)) {
          std::cout << "group1" << endl;

          triangulate_singularity_case2(singular_point, singular_direction,
                                        &my_mesh,
                                        (i == 0) ? kInvalidPointIndex : 0);
        } else if (std::regex_match(name, A2) ) {
          std::cout << "group2" << endl;
          triangulate_singularity_circular(singular_point, singular_direction,
                                           &my_mesh,
                                           (i == 0) ? kInvalidPointIndex : 0);
        } else if (std::regex_match(name, A3)) {
          std::cout << "group3" << endl;
          triangulate_An_analytical(singular_point, singular_direction,
                                    &my_mesh, (i == 0) ? kInvalidPointIndex : 0,
                                    types[j]);
        } else {
          std::cout << "group4" << endl;
          triangulate_cone_iterative(singular_point, singular_direction,
                                     &my_mesh,
                                     (i == 0) ? kInvalidPointIndex : 0);
        }
      }
    }
  }
*/
  explicit BasicAlgorithm(string name, Function f, const Function &F,
                          const Function &G, const vector<Point> &polyline,
                          numeric e_size, realsymbol x, realsymbol y,
                          realsymbol z, BoundingBox bounding_box)
      : name(name),
        F(f),
        e_size(e_size),
        x(x),
        y(y),
        z(z),
        bounding_box(bounding_box) {
    get_local_mesh_curve(&my_mesh, F, G, f, polyline, bounding_box, e_size);
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
  int fix_holes(const HalfEdge &working_edge);
  bool fix_holes2(const HalfEdge &working_edge);
  bool fix_breakers(const HalfEdge &working_edge, const Point &projected,
                    const bool Delaunay);
  bool fix_same_points(const HalfEdge &working_edge, const Point &projected);

  pair<MeshPoint, MeshPoint> find_prev_next(const HalfEdge &working_edge) const;
  MeshPointIndex find_triangle_hole(const HalfEdge &working_edge);
  std::pair<FaceIndex, FaceIndex> find_neighbor_faces(
      const HalfEdge &working_edge, const Point &P) const;
  bool neighbor_triangles_normal_check(const Triangle &T1,
                                       const Triangle &T2) const;
  bool orientability_check(const HalfEdge &working_edge,
                           const MeshPoint &point) const;

  std::optional<Point> get_projected(const HalfEdge &working_edge) const;
  bool overlap_normals_check(const MeshPoint &candidate,
                             const HalfEdge &working_edge) const;
  std::optional<vector<MeshPoint>> find_close_points(
      const Point &P, const HalfEdge &working_edge) const;
  int update_border(const Edge &new_edge1, const Edge &new_edge2);
  bool check_conditions(const HalfEdge &working_edge, const Point &P,
                        const bool Delaunay,
                        const MeshPointIndex i_P = kInvalidPointIndex) const;
  bool Delaunay_conditions(const HalfEdge &working_edge, const Point &P) const;
  bool Delaunay_conditions_debug(const HalfEdge &working_edge,
                                 const Point &P) const;
  bool non_Delaunay_conditions(const HalfEdge &working_edge, const Point &P,
                               const MeshPointIndex i_P) const;

  void create_triangle(const HalfEdge &working_edge, const Point &P,
                       const bool is_new, const MeshPointIndex index_P = -1);
  bool gc_distance_check(const Triangle T) const;
  bool good_edges(const HalfEdge &working_edge, const Point &P) const;
  bool good_new_point(const MeshPoint &point) const;
  bool basic_triangle(const HalfEdge &working_edge, const MeshPoint &point);

  bool is_active(const HalfEdgeIndex &halfedge_index) const;
  bool is_checked(const HalfEdgeIndex &halfedge_index) const;
  bool is_bounding(const HalfEdgeIndex &halfedge_index) const;
  bool is_border(const HalfEdgeIndex &halfedge_index) const;
  // bool is_border_point(const MeshPoint &P) const;
  bool is_border_point(const Point &P) const;
  void delete_from_active(const HalfEdgeIndex &halfedge_index);
  void delete_from_checked(const HalfEdgeIndex &halfedge_index);
  void push_edge_to_active(const HalfEdgeIndex &halfedge_index);
  void push_edge_to_checked(const HalfEdgeIndex &halfedge_index);

  std::optional<MeshPoint> get_closest_point(
      const HalfEdge &working_edge) const;
  std::optional<pair<Edge, numeric>> get_closest_edge(const Point &P,
                                                      const Triangle &N) const;

  void fix_corners();

  // section seed triangle
  Edge get_seed_edge(Point seed_point) const;
  Point get_seed_triangle(const Edge &e) const;
  Triangle find_seed_triangle(Point seed) const;
  void create_mesh(Point seed);

  // section singularities
  void triangulate_A1_starter(const Point &singular,
                              const Vector &singular_direction, Mesh *mesh,
                              const MeshPointIndex singular_index);
  void triangulate_singularity_circular(const Singularity &singularity,
                                        const int branch,
                                        const MeshPointIndex singular_index);
  void triangulate_An_analytical(const Singularity &singularity,
                                 const int branch,
                                 const MeshPointIndex singular_index);
  void triangulate_singularity_case2(const Singularity &singularity,
                                     const int branch,
                                     const MeshPointIndex singular_index);
  void triangulate_cone_iterative(const Point &singular,
                                  const Vector &singular_direction, Mesh *mesh,
                                  const MeshPointIndex singular_index);

  void get_local_mesh_curve(Mesh *local_mesh, const Function &_F,
                            const Function &_G, const Function &inter_FG,
                            const vector<Point> &polyline,
                            const BoundingBox &bounding_box,
                            const numeric e_size);
  void get_local_mesh_point(const vector<vector<Point>> &points,
                            const Singularity &singularity,
                            const MeshPointIndex singular_index,
                            const int num_triangles, const int branch,
                            const int layers);
  void get_points_singular_point(const Singularity singularity,
                                 const int num_triangles,
                                 const numeric &edge_length,
                                 vector<vector<Point>> &points,
                                 const int branch = 0, const int layers = 1);

 private:
  string name;
  const Function F;
  // vector<Edge> active_edges;
  // vector<Edge> checked_edges;
  // vector<Edge> bounding_edges;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
  BoundingBox bounding_box;
  HalfEdgeIndex _find_prev(const HalfEdge &working_edge) const;
  HalfEdgeIndex _find_next(const HalfEdge &working_edge) const;
  bool _fix_prev_next(const HalfEdge &working_edge, const MeshPoint &prev,
                      const bool Delaunay);
  std::optional<HalfEdgeIndex> _find_closest_prev(
      const HalfEdge &working_edge, const vector<HalfEdgeIndex> &prev) const;
  std::optional<HalfEdgeIndex> _find_closest_next(
      const HalfEdge &working_edge, const vector<HalfEdgeIndex> &next) const;
  std::optional<vector<MeshPoint>> _linear_close_points_finder(
      const Point &P, const HalfEdge &working_edge) const;
  std::optional<vector<MeshPoint>> _tree_close_points_finder(
      const Point &P, const HalfEdge &working_edge) const;
  bool _check_normals(const HalfEdge &working_edge, const Point &point) const;
};
