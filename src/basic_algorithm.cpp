#include "basic_algorithm.h"

#include <iostream>

using std::cout;

/*
void BasicAlgorithm::add_marks() {
  auto border = connect_edges(active_edges, checked_edges);
  for (auto edge : border) {
    auto dir = e_size / 5 *
               find_direction(edge, my_mesh.find_triangle_with_edge(edge));
    assertm(!Vector(edge.A(), edge.B()).is_zero(), "Zero edge vector!");
    Vector edge_dir = e_size / 5 * Vector(edge.A(), edge.B()).unit();
    Edge new_e(Point(edge.get_midpoint(), edge_dir.vector_inverse() / 2),
               Point(edge.get_midpoint(), edge_dir / 2));
    auto new_p = Point(edge.get_midpoint(), dir);
    my_mesh.add_triangle(new_e, new_p);
  }
}
*/
// returns vector of points closer than 0.4*e_size to point P sorted from
// closest to working edge if there are no points returns std::nullopt
std::optional<vector<MeshPoint>> BasicAlgorithm::find_close_points(
    const Point &P, const HalfEdge &working_edge,
    const Face &incident_face) const {
  numeric min_dist = 0.4 * e_size;
  vector<MeshPoint> close_points;
  for (MeshPoint meshpoint : my_mesh.get_mesh_points()) {
    if (my_mesh.is_boundary_point(meshpoint) &&
        Vector(meshpoint.get_point(), P).get_length() < min_dist)
      close_points.push_back(meshpoint);
  }

  /*
    vector<Edge> edges =
        connect_edges(connect_edges(my_mesh.get_active_edges(),
    my_mesh.get_checked_edges()), my_mesh.get_boundary_edges()); for (auto edge
    : edges) { if (edge.A() != P && edge.A() != working_edge.A() && edge.A() !=
    working_edge.B()) { numeric dist = Vector(edge.A(), P).get_length(); if
    (dist < min_dist) { bool found = false; for (auto point : close_points) { if
    (point == edge.A()) { found = true;
            }
          }
          if (!found) {
            close_points.push_back(edge.A());
          }
        }
      }
      if (edge.B() != P && edge.B() != working_edge.A() &&
          edge.B() != working_edge.B()) {
        numeric dist = Vector(edge.B(), P).get_length();
        if (dist < min_dist) {
          bool found = false;
          for (auto point : close_points) {
            if (point == edge.B()) {
              found = true;
            }
          }
          if (!found) {
            close_points.push_back(edge.B());
          }
        }
      }
    }
    */
  std::cout << "close points count: " << close_points.size() << endl;
  if (close_points.empty()) return std::nullopt;

  sort(close_points.begin(), close_points.end(),
       [this, &working_edge, &incident_face](auto i, auto j) {
         return line_point_dist(my_mesh, working_edge, i, incident_face) <
                line_point_dist(my_mesh, working_edge, j, incident_face);
       });

  return close_points;
}
/*
// fixes corners of bounded triangulation
void BasicAlgorithm::fix_corners() {
  realsymbol my_x("my_x"), my_y("my_y"), my_z("my_z");
  for (auto edge : bounding_edges) {
    // binary number, ones at place of side at which point is lying on
    int faces_values_A = bounding_box.faces(edge.A());
    int faces_values_B = bounding_box.faces(edge.B());
    // logical and
    int common_faces = faces_values_A & faces_values_B;
    int index_A;
    int index_B;
    if (common_faces != 0) continue;
    for (int i = 0; i < 6; ++i) {
      // value of ith bit
      if (faces_values_A & (1 << i)) {
        index_A = i;
      }
      if (faces_values_B & (1 << i)) {
        index_B = i;
      }
    }
    // direciton vector of intersection line of two faces on which the edge lies
    Vector v(1, 1, 1);
    std::optional<numeric> Px;
    std::optional<numeric> Py;
    std::optional<numeric> Pz;
    if (index_A == 0 || index_B == 0 || index_A == 1 || index_B == 1) {
      v = v - Vector(1, 0, 0);
      if (index_A == 0 || index_B == 0) {
        Px = bounding_box.min_x();
      } else if (index_A == 1 || index_B == 1) {
        Px = bounding_box.max_x();
      } else {
        assertm(false, "Point of an edge lying on x_min and x_max sides!");
      }
    }
    if (index_A == 2 || index_B == 2 || index_A == 3 || index_B == 3) {
      v = v - Vector(0, 1, 0);
      if (index_A == 2 || index_B == 2) {
        Py = bounding_box.min_y();
      } else if (index_A == 3 || index_B == 3) {
        Py = bounding_box.max_y();
      } else {
        assertm(false, "Point of an edge lying on y_min and y_max sides!");
      }
    }
    if (index_A == 4 || index_B == 4 || index_A == 5 || index_B == 5) {
      v = v - Vector(0, 0, 1);
      if (index_A == 4 || index_B == 4) {
        Pz = bounding_box.min_z();
      } else if (index_A == 5 || index_B == 5) {
        Pz = bounding_box.max_z();
      } else {
        assertm(false, "Point of an edge lying on z_min and z_max sides!");
      }
    }
    if (!Px.has_value()) Px = edge.get_midpoint().x();
    if (!Py.has_value()) Py = edge.get_midpoint().y();
    if (!Pz.has_value()) Pz = edge.get_midpoint().z();

    std::optional<Point> projected = std::nullopt;
    Point P(Px.value(), Py.value(), Pz.value());
    if (v.is_zero()) {
      projected = P;
    } else {
      projected = project(P, v, F, e_size);
    }
    assertm(projected.has_value(), "Point without value!");
    if (Vector(projected.value(), edge.get_midpoint()).get_length() <
        edge.get_length()) {
      if(Triangle(edge.A(), edge.B(), projected.value()).is_triangle()){ // &&
good_orientation(edge, projected, my_mesh.find_triangle_with_edge(edge))
        create_triangle(edge, projected.value());
      }
    }
  }
  return;
}
*/

// checks if conditions required in the first part of the algorithm are
// satisfied
bool BasicAlgorithm::Delaunay_conditions(const HalfEdge &working_edge,
                                         const Point &P,
                                         const Face &incident_face) const {
  if (working_edge.get_point_A() == P || working_edge.get_point_B() == P)
    return false;
  Triangle T =
      Triangle(working_edge.get_point_B(), working_edge.get_point_A(), P);
  // auto [prev, next] = find_prev_next(working_edge, incident_face);
  // all the conditions for a new triangle in the first part of algorithm
  bool is_triangle = T.is_triangle();
  bool has_good_orientation =
      good_orientation(my_mesh, working_edge, P, incident_face);
  bool delaunay =
      my_mesh.check_Delaunay(my_mesh, T, working_edge, incident_face);
  bool has_good_edges = good_edges(working_edge, P);
  if (!is_triangle) cout << "Triangle condition failed!";
  if (!has_good_orientation) cout << "Orientation condition failed!";
  if (!delaunay) cout << "Delaunay condition failed!";
  if (!has_good_edges) cout << "Good Edges condition failed!";

  return (is_triangle && delaunay /* && has_good_edges*/);
}

// checks if conditions required in the first part of the algorithm are
// satisfied
bool BasicAlgorithm::non_Delaunay_conditions(const HalfEdge &working_edge,
                                             const Point &P,
                                             const Face &incident_face) const {
  if (working_edge.get_point_A() == P || working_edge.get_point_B() == P)
    return false;
  Triangle T =
      Triangle(working_edge.get_point_A(), working_edge.get_point_B(), P);
  auto [prev_point, next_point] = find_prev_next(working_edge, incident_face);
  // all the conditions for a new triangle in the second part of algorithm

  bool is_triangle = T.is_triangle();
  bool has_good_orientation =
      good_orientation(my_mesh, working_edge, P, incident_face);
  bool has_good_edges = good_edges(working_edge, P);
  if (!is_triangle) cout << "Triangle condition failed!";
  if (!has_good_orientation) cout << "Orientation condition failed!";
  if (!has_good_edges) cout << "Good Edges condition failed!";

  return (is_triangle /*&& has_good_edges*/);
}

// creates new triangle and adds it to mesh
void BasicAlgorithm::create_triangle(const HalfEdge &working_edge,
                                     const Point &P, const std::string type,
                                     const MeshPointIndex index_P) {
  // std::cout << "in create_triangle: " << endl
  //           << "working_edge index is: " << working_edge.get_index() << endl;
  //  my_mesh.cout_mesh();
  assertm(working_edge.get_point_A() != P && working_edge.get_point_B() != P,
          "Wrong new triangle!");
  Triangle new_triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                        P);
  assertm(new_triangle.is_triangle(), "New triangle is not a triangle!");
  my_mesh.add_triangle(working_edge.get_index(), P, type, index_P);
  // my_mesh.cout_mesh();
  //  update_border(new_edge1, new_edge2);
  return;
}

// updates active and checked edges and returns number of new edges
/*
int BasicAlgorithm::update_border(const Edge &new_edge1,
                                  const Edge &new_edge2) {
  assertm(new_edge1 != new_edge2, "Same edges while updating border!");
  int new_edges = 0;

  if (is_border(new_edge1)) {
    delete_from_active(new_edge1);
    delete_from_checked(new_edge1);
  } else {
    push_edge_to_active(new_edge1);
    new_edges++;
  }

  if (is_border(new_edge2)) {
    delete_from_active(new_edge2);
    delete_from_checked(new_edge2);
  } else {
    push_edge_to_active(new_edge2);
    new_edges++;
  }

  if (bounding_box.new_bounding_edge(new_edge1)) {
    delete_from_active(new_edge1);
    delete_from_checked(new_edge1);
    if (!is_bounding(new_edge1)) bounding_edges.push_back(new_edge1);
  }
  if (bounding_box.new_bounding_edge(new_edge2)) {
    delete_from_active(new_edge2);
    delete_from_checked(new_edge2);
    if (!is_bounding(new_edge2)) bounding_edges.push_back(new_edge2);
  }

  return new_edges;
}
*/

// makes triangle if prev == next
// TODO basic triangle na okraji!!
bool BasicAlgorithm::basic_triangle(const HalfEdge &working_edge,
                                    const Face &incident_face,
                                    const MeshPoint &prev,
                                    const MeshPoint &next) {
  // determine other point on the neighbour triangle - opposite to working_edge
  MeshPoint opposite_point = my_mesh.get_meshpoint(
      my_mesh.get_halfedge(working_edge.get_next()).get_B());

  // if prev and next are the same and its not only one triangle create triangle
  if (prev == next && opposite_point != prev &&
      Triangle(working_edge.get_point_B(), working_edge.get_point_A(),
               prev.get_point())
          .is_triangle()) {
    create_triangle(working_edge, prev, "fill", prev.get_index());
    return true;
  }
  return false;
}

// true if halfedge is active
bool BasicAlgorithm::is_active(const HalfEdgeIndex &halfedge_index) const {
  return my_mesh.is_active(halfedge_index);
}
// true if halfedge is checked
bool BasicAlgorithm::is_checked(const HalfEdgeIndex &halfedge_index) const {
  return my_mesh.is_checked(halfedge_index);
}
// true if halfedge is active
bool BasicAlgorithm::is_bounding(const HalfEdgeIndex &halfedge_index) const {
  return my_mesh.is_bounding(halfedge_index);
}
// true if halfedge is active or checked
bool BasicAlgorithm::is_border(const HalfEdgeIndex &halfedge_index) const {
  return (is_active(halfedge_index) || is_checked(halfedge_index) ||
          is_bounding(halfedge_index));
}

// true if point is on border of mesh
bool BasicAlgorithm::is_border_point(const MeshPoint &P) const {
  return my_mesh.is_boundary_point(P);
}

bool BasicAlgorithm::is_border_point(const Point &P) const {
  for (auto point : my_mesh.get_mesh_points()) {
    if (is_border_point(point) &&
        Vector(point.get_point(), P).get_length() < kEps)
      return true;
  }
}

// throws error if it is found more than once
void BasicAlgorithm::delete_from_active(const HalfEdgeIndex &halfedge_index) {
  assertm(is_active(halfedge_index), "Deleting non-active edge!");
  my_mesh.remove_active_edge(halfedge_index);
}

// throws error if it is found more than once
void BasicAlgorithm::delete_from_checked(const HalfEdgeIndex &halfedge_index) {
  assertm(is_checked(halfedge_index), "Deleting non-checked edge!");
  my_mesh.remove_checked_edge(halfedge_index);
}

// throws error if it is already there
void BasicAlgorithm::push_edge_to_active(const HalfEdgeIndex &halfedge_index) {
  assertm(!is_active(halfedge_index), "Edge already in active edges!");
  my_mesh.add_edge_to_active(halfedge_index);
  return;
}

// throws error if it is already there
void BasicAlgorithm::push_edge_to_checked(const HalfEdgeIndex &halfedge_index) {
  assertm(!is_checked(halfedge_index), "Edge already in checked_edges!");
  my_mesh.add_edge_to_checked(halfedge_index);
  return;
}

// finds closest border point to edge
std::optional<Point> BasicAlgorithm::get_closest_point(
    const HalfEdge &working_edge, const Face &incident_face) const {
  std::vector<MeshPoint> meshpoints = my_mesh.get_mesh_points();
  std::optional<MeshPoint> closest_point = std::nullopt;
  std::optional<numeric> min_dist;

  for (auto meshpoint : meshpoints) {
    auto dist = line_point_dist(my_mesh, working_edge, meshpoint.get_point(),
                                incident_face);
    if (dist < min_dist &&
        meshpoint.get_point() != working_edge.get_point_A() &&
        meshpoint.get_point() != working_edge.get_point_B() &&
        non_Delaunay_conditions(working_edge, meshpoint.get_point(),
                                incident_face))
      min_dist = dist;
    closest_point = meshpoint;
  }
  if (closest_point.has_value() && min_dist.has_value() &&
      min_dist.value() < 1.5 * e_size) {
    return closest_point.value();
  }
  return std::nullopt;

  /*
  auto border_edges = connect_edges(active_edges, checked_edges);
  std::optional<Point> closest_point = std::nullopt;
  std::optional<numeric> min_dist;
  for (auto edge : border_edges) {
    if (!closest_point.has_value()) {
      if (edge.A() != working_edge.A() && edge.A() != working_edge.B() &&
          non_Delaunay_conditions(working_edge, edge.A(), neighbour_triangle)) {
        closest_point = edge.A();
        min_dist = line_point_dist(working_edge, closest_point.value(),
                                   neighbour_triangle);
      }
      if (!closest_point.has_value() && edge.B() != working_edge.A() &&
          edge.B() != working_edge.B() &&
          non_Delaunay_conditions(working_edge, edge.B(), neighbour_triangle)) {
        closest_point = edge.B();
        min_dist = line_point_dist(working_edge, closest_point.value(),
                                   neighbour_triangle);
      }
    }
    if (closest_point.has_value()) {
      if (line_point_dist(working_edge, edge.A(), neighbour_triangle) <
              min_dist &&
          edge.A() != working_edge.A() && edge.A() != working_edge.B() &&
          non_Delaunay_conditions(working_edge, edge.A(), neighbour_triangle)) {
        closest_point = edge.A();
        min_dist = line_point_dist(working_edge, closest_point.value(),
                                   neighbour_triangle);
      }
      if (line_point_dist(working_edge, edge.B(), neighbour_triangle) <
              min_dist &&
          edge.B() != working_edge.A() && edge.B() != working_edge.B() &&
          non_Delaunay_conditions(working_edge, edge.B(), neighbour_triangle)) {
        closest_point = edge.B();
        min_dist = line_point_dist(working_edge, closest_point.value(),
                                   neighbour_triangle);
      }
    }
  }

  if (closest_point.has_value() && min_dist.has_value() &&
      min_dist.value() < 1.5 * e_size) {
    return closest_point.value();
  }
  return std::nullopt;
  */
}

// checks if edges of new triangle are active or are not im mesh
bool BasicAlgorithm::good_edges(const HalfEdge &working_edge,
                                const Point &P) const {
  Edge new_edge1(P, working_edge.get_point_A());
  Edge new_edge2(working_edge.get_point_B(), P);

  HalfEdgeIndex new_edge1_index = my_mesh.get_edge_index(new_edge1);
  HalfEdgeIndex new_edge2_index = my_mesh.get_edge_index(new_edge2);

  // Edges are good if either they are not in mesh
  // or they are in mesh and are border

  bool is_in_mesh_1 = (new_edge1_index != kInvalidEdgeIndex);
  bool is_in_mesh_2 = (new_edge2_index != kInvalidEdgeIndex);

  if (!is_in_mesh_1 && !is_in_mesh_2) return true;

  bool is_good_1 = true, is_good_2 = true;

  if (is_in_mesh_1) is_good_1 = is_border(new_edge1_index);
  if (is_in_mesh_2) is_good_2 = is_border(new_edge2_index);

  return (is_good_1 && is_good_2);
}
/*
// finds closest border edge to point P
std::optional<pair<Edge, numeric>> BasicAlgorithm::get_closest_edge(
    const Point &P, const Triangle &N) const {
  auto border = connect_edges(active_edges, checked_edges);
  std::optional<pair<Edge, numeric>> closest_edge = std::nullopt;
  numeric dist = 0;
  for (auto edge : border) {
    dist = line_point_dist(edge, P, N);
    if (!closest_edge.has_value()) {
      closest_edge = pair(edge, dist);
    } else if (dist < closest_edge.value().second) {
      closest_edge = pair(edge, dist);
    }
  }
  assertm(closest_edge.has_value(), "Edge without value!");
  return closest_edge;  // TODO: return non-optional
}
*/

bool BasicAlgorithm::overlap_normals_check(const MeshPoint &candidate,
                                           const HalfEdge &working_edge) const {
  // if (candidate == prev || candidate == next || candidate == working_edge.A()
  // ||
  //     candidate == working_edge.B())
  //   return false;

  if (candidate.get_point() == working_edge.get_point_A() ||
      candidate.get_point() == working_edge.get_point_B())
    return false;

  Triangle my_triangle(working_edge.get_point_A(), working_edge.get_point_B(),
                       candidate.get_point());

  if (my_triangle.is_triangle() /*&&
      good_edges(working_edge, candidate.get_point())*/) {
    Vector my_normal = F.outside_normal(my_triangle, e_size);
    Face overlap_face = my_mesh.get_face(
        my_mesh.get_halfedge(candidate.get_outgoing()[0]).get_incident());
    Vector overlap_normal =
        F.outside_normal(overlap_face.get_triangle(), e_size);
    // if normals have the same orientation
    if (overlap_normal * my_normal > 0) {
      return true;
    }
    return false;
  }
  return false;
}

Point BasicAlgorithm::get_projected(const HalfEdge &working_edge,
                                    const Face &incident_face) const {
  // cout << "in get_projected" << endl;
  Point center = working_edge.get_midpoint();
  assertm(Vector(working_edge.get_point_A(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6 &&
              Vector(working_edge.get_point_B(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6,
          "Wrong get_midpoint function!");

  auto [prev_point, next_point] = find_prev_next(working_edge, incident_face);

  // height of equilateral triangle based on e_size
  numeric basic_height = e_size * sqrt(numeric(3)) / 2;

  numeric average_edge_length =
      (1 / numeric(3)) *
      (Edge(working_edge.get_point_A(), prev_point.get_point()).get_length() +
       Edge(working_edge.get_point_B(), next_point.get_point()).get_length() +
       working_edge.get_length());

  // non adaptive height
  numeric height = basic_height;

  // height of equilateral triangle based on working_edge size
  // numeric height = working_edge.get_length() * sqrt(numeric(3)) / 2;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;

  // height of equilateral triangle based on neighbour edges size
  // numeric height = average_edge_length * sqrt(numeric(3)) / 2;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;

  // height of equilateral triangle based on neighbour edges size with influence
  // of e_size numeric height = 0.75*(average_edge_length * sqrt(numeric(3)) /
  // 2) + 0.25*basic_height; if(height<e_size/3) height = 0.25*(e_size/3) +
  // 0.75*height; if(height>2*e_size) height = 2*e_size;

  Vector direction = height * find_direction(working_edge, incident_face);
  assertm(direction * incident_face.get_normal() < kEps, "Wrong direction!");
  assertm(direction * Vector(working_edge.get_point_A(),
                             working_edge.get_point_B()) <
              kEps,
          "Wrong direction!");

  Point P(center, direction);
  assertm(Vector(center, P).get_length() - height < kEps,
          "Wrong point to project!");

  /*
  Vector n_A = F.get_gradient_at_point(working_edge.A()).unit();
  Vector n_B = F.get_gradient_at_point(working_edge.B()).unit();

  Vector normal = (n_A + n_B) / 2;
  */

  assertm(!F.get_gradient_at_point(P).is_zero(), "Zero gradient!");
  Vector normal = F.get_gradient_at_point(P).unit();
  Point projected = project(P, normal, F, {e_size});

  /*
  assertm((Vector(working_edge.A(), projected).get_length() < 4 * e_size) &&
              (Vector(working_edge.B(), projected).get_length() < 4 * e_size),
          "Projected point too far!");
  */
  return projected;
}

pair<MeshPoint, MeshPoint> BasicAlgorithm::find_prev_next(
    const HalfEdge &working_edge, const Face &incident_face) const {
  auto [prev_edge_index, next_edge_index] = pair<HalfEdgeIndex, HalfEdgeIndex>(
      _find_prev(working_edge, incident_face),
      _find_next(working_edge, incident_face));
  const HalfEdge prev_edge = my_mesh.get_halfedge(prev_edge_index);
  const HalfEdge next_edge = my_mesh.get_halfedge(next_edge_index);
  const MeshPoint prev_point = my_mesh.get_meshpoint(prev_edge.get_A());
  const MeshPoint next_point = my_mesh.get_meshpoint(next_edge.get_B());
  return {prev_point, next_point};
}

// returns true if prev/next triangle is added to mesh else returns false
bool BasicAlgorithm::fix_prev_next(const HalfEdge &working_edge,
                                   const Face &incident_face,
                                   const bool is_prev, const bool Delaunay) {
  // find previous and next points
  auto [prev_point, next_point] = find_prev_next(working_edge, incident_face);

  // assertm(Vector(working_edge.A(), prev).get_length() < 3 * e_size,
  //         "Wrong prev point!");
  // assertm(Vector(working_edge.B(), next).get_length() < 3 * e_size,
  //         "Wrong prev point!");

  // assertm(is_border(Edge(working_edge.get_point_A(), prev)), "Prev edge not
  // in
  //  border!");
  //  assertm(is_border(Edge(working_edge.B(), next)), "Next edge not in
  //  border!");
  // cout << "in fixprevnext" << endl;
  MeshPoint vertex = prev_point;
  Edge edge = working_edge.get_edge();
  // fix so that vertex is the new vertex and edge.A() is the adjacent vertex
  if (!is_prev) {
    edge = Edge(working_edge.get_point_B(), working_edge.get_point_A());
    vertex = next_point;
  }

  // potentialy new triangle
  Triangle potential_triangle(edge.A(), edge.B(), vertex.get_point());

  if (potential_triangle.AB().get_length() > 2 * e_size ||
      potential_triangle.CA().get_length() > 2 * e_size ||
      potential_triangle.BC().get_length() > 2 * e_size) {
    return false;
  }

  // checks if the potential triangle has good orientation and angle near A is
  // less than 90 degrees and checks Delaunay
  bool del_con =
      Delaunay_conditions(working_edge, vertex.get_point(), incident_face);
  if ((/*Delaunay **/  del_con 
  /*|| !Delaunay * non_Delaunay_conditions(working_edge, vertex.get_point(),
                                           incident_face)*/)) {
    assertm(prev_point != next_point, "Prev and next are the same!");
    // cout << "creating prev triangle!" << endl;
    create_triangle(working_edge, vertex.get_point(),
                    is_prev ? "previous" : "next", vertex.get_index());
    return true;
  }
  return false;
}

// returns true if overlap triangle is added to mesh else returns false
bool BasicAlgorithm::fix_overlap(const HalfEdge &working_edge,
                                 const Face &incident_face,
                                 const MeshPoint &overlap_point,
                                 const bool Delaunay) {
  auto [prev_point, next_point] = find_prev_next(working_edge, incident_face);

  // cout << "prev: " << prev_point << endl;
  // cout << "next: " << next_point << endl;
  // cout << "overlap: " << overlap_point << endl;
  if (overlap_point == prev_point && overlap_point == next_point) {
    // cout << "basic" << endl;
    return basic_triangle(working_edge, incident_face, prev_point, next_point);
  }
  if (overlap_point == prev_point) {
    // cout << "prev" << endl;
    return fix_prev_next(working_edge, incident_face, true, true);
  }
  if (overlap_point == next_point) {
    // cout << "next" << endl;
    return fix_prev_next(working_edge, incident_face, false, true);
  }

  // checks if overlap point is not neighbour or working edge point also if
  // overlap is on the border and if overlap triangle has good orientation of
  // normal
  if (overlap_normals_check(overlap_point, working_edge)) {
    Triangle potential_triangle =
        Triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                 overlap_point.get_point());

    // if Delaunay constraint is satisfied add the triangle to
    // triangulation and end
    if (Delaunay * Delaunay_conditions(working_edge, overlap_point.get_point(),
                                       incident_face) ||
        !Delaunay * non_Delaunay_conditions(working_edge,
                                            overlap_point.get_point(),
                                            incident_face)) {
      // assertm(Vector(working_edge.A(), overlap_point).get_length() <
      //             3 * working_edge.get_length(),
      //         "Weird distance of overlap point!");
      // assertm(Vector(working_edge.B(), overlap_point).get_length() <
      //             3 * working_edge.get_length(),
      //         "Weird distance of overlap point!");

      assertm(overlap_point.get_point() != prev_point.get_point() &&
                  overlap_point.get_point() != next_point.get_point(),
              "Creating overlap from non overlap point!");
      create_triangle(working_edge, overlap_point.get_point(), "overlap",
                      overlap_point.get_index());

      // cout << "overlap" << endl;
      return true;
    }
  }
  return false;
}

bool BasicAlgorithm::fix_proj(const HalfEdge &working_edge,
                              const Point &projected, const Face &incident_face,
                              const MeshPoint &prev, const MeshPoint &next) {
  Triangle potential_triangle(working_edge.get_point_A(),
                              working_edge.get_point_B(), projected);
  assertm(potential_triangle.is_triangle(), "Projected triangle not valid!");

  // checks if there are some points very close to projected point
  // as surrounding points were taken working_edge points

  if (auto surrounding_points =
          find_close_points(projected, working_edge, incident_face);
      surrounding_points.has_value()) {
    // points closer to projected point than 0.4*e_size sorted from closest
    vector<MeshPoint> close_points = surrounding_points.value();
    for (auto close_point : close_points) {
      if (Delaunay_conditions(working_edge, close_point.get_point(),
                              incident_face)) {
        // if close point is prev we want to try fix prev
        if (close_point == prev) {
          if (fix_prev_next(working_edge, incident_face, true, true)) {
            return true;
          }
        }
        // if close point is next we want to try fix next
        else if (close_point == next) {
          if (fix_prev_next(working_edge, incident_face, false, true)) {
            return true;
          }
        }

        // if close point is overlap we want to try fix overlap
        else {
          if (fix_overlap(working_edge, incident_face, close_point, true)) {
            return true;
          }
        }
      }
    }
    // if there are close points but nothing worked we want to try to
    // construct original triangle
  }

  /*
  auto close_edge = get_closest_edge(projected, neighbour_triangle).value();
  assertm(is_border(close_edge.first), "Closest edge not in border edges!");
  if (close_edge.second < e_size / 3) {
    Edge closest_edge = close_edge.first;
    Point P1 = closest_edge.A();
    Point P2 = closest_edge.B();


    //Vector n_A = F.get_gradient_at_point(closest_edge.A()).unit();
    //Vector n_B = F.get_gradient_at_point(closest_edge.B()).unit();
    //Vector normal = (n_A + n_B) / 2;


    assertm(!F.get_gradient_at_point(closest_edge.get_midpoint()).is_zero(),
            "Zero gradient!");
    Vector normal = F.get_gradient_at_point(closest_edge.get_midpoint()).unit();

    if (P1 != working_edge.A() && P1 != working_edge.B()) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P1);
      if (Delaunay_conditions(working_edge, P1, neighbour_triangle)) {
        create_triangle(working_edge, P1);
        { return true; }
      }
    }

    if (P2 != working_edge.A() && P2 != working_edge.B()) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P2);
      if (Delaunay_conditions(working_edge, P2, neighbour_triangle)) {
        create_triangle(working_edge, P2);
        { return true; }
      }
    }
    // TODO opraviť divide_triangle_by_point

    // point in the middle of side
    Point P3 = project(closest_edge.get_midpoint(), normal, F, {e_size});
    // Point P3 = closest_edge.get_midpoint();
    if (P3 != working_edge.A() && P3 != working_edge.B()) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P3);
      if (Delaunay_conditions(working_edge, P3, neighbour_triangle)) {
        // point P3 is not a vertex, so we need to subdivide triagle with this
        // edge and add subdivided triangles to mesh

        Edge new_edge1(working_edge.A(), P3);
        Edge new_edge2(working_edge.B(), P3);

        delete_from_active(closest_edge);
        delete_from_checked(closest_edge);
        my_mesh.divide_triangle_by_point(closest_edge, P3);
        // update_border(new_edge1, new_edge2);
        // if(working_edge.A() != closest_edge.A() && working_edge.B() !=
        // closest_edge.A())
        push_edge_to_active(Edge(closest_edge.A(), P3));
        // if(working_edge.A() != closest_edge.B() && working_edge.B() !=
        // closest_edge.B())
        push_edge_to_active(Edge(closest_edge.B(), P3));

        my_mesh.add_triangle(working_edge, P3);

        update_border(new_edge1, new_edge2);
        return true;
      }
    }
  }
  */
  if (Delaunay_conditions(working_edge, projected, incident_face)
   /* my_mesh.check_Delaunay(maybe_new_T) &&
      good_edges(working_edge, projected)*/) {
    /*Point clipped = bounding_box.crop_to_box(working_edge.get_midpoint(),
                                             projected, e_size, F);*/
    // if (non_Delaunay_conditions(working_edge, clipped, neighbour_triangle)) {
    // if (non_Delaunay_conditions(working_edge, projected, incident_face)) {
    // create_triangle(working_edge, clipped);
    assertm(projected != prev.get_point() && projected != next.get_point(),
            "Creating new triangle with prev or next!");
    create_triangle(working_edge, projected, "new");
    return true;
    //}
  }

  // if nothing worked
  return false;
}

bool BasicAlgorithm::fix_breakers(const HalfEdge &working_edge,
                                  const Point &projected,
                                  const Face &incident_face,
                                  const bool Delaunay) {
  Triangle proj_T(working_edge.get_point_A(), working_edge.get_point_B(),
                  projected);

  // points that break Delaunay constraint
  vector<MeshPoint> breakers = my_mesh.get_breakers(proj_T);
  std::cout << "breakers count: " << breakers.size() << endl;
  // sort from closest to working_edge
  std::sort(breakers.begin(), breakers.end(),
            [this, &working_edge, &incident_face](auto i, auto j) {
              return line_point_dist(my_mesh, working_edge, i, incident_face) <
                     line_point_dist(my_mesh, working_edge, j, incident_face);
            });

  // try create triangle with breakers

  for (auto point : breakers) {
    if (is_border_point(point) &&
        fix_overlap(working_edge, incident_face, point, Delaunay)) {
      // cout << "here" << endl;
      return true;
    }
  }

  return false;
}

// one step of the algorithm
bool BasicAlgorithm::step(const HalfEdgeIndex &working_edge_index) {
  const FaceIndex incident_face_index =
      my_mesh.get_incident_face(working_edge_index);
  const HalfEdge working_edge = my_mesh.get_halfedge(working_edge_index);
  const Face incident_face = my_mesh.get_face(incident_face_index);

  auto [prev_point, next_point] = find_prev_next(working_edge, incident_face);
  // find candidate point for working_edge
  Point projected = get_projected(working_edge, incident_face);

  // cout << "Incident face: " << endl << incident_face.get_triangle() << endl;
  // cout << "Projected: " << projected << endl;

  assertm(!is_active(working_edge_index),
          "Working edge found in active edges!");
  assertm(!is_checked(working_edge_index),
          "Working edge found in checked edges!");
  assertm(!is_bounding(working_edge_index),
          "Working edge found in bounding edges!");
  assertm(!is_border(working_edge_index),
          "Working edge found in border edges!");

  /*
    // if point projects on already existing border point
    if (fix_same_points(working_edge, projected, neighbour_triangle)) {
      return true;
    }
*/
  // if there is a hole in triangle shape, fill it with triangle
  if (basic_triangle(working_edge, incident_face, prev_point, next_point)) {
    cout << "1" << endl;
    return true;
  }

  if (fix_breakers(working_edge, projected, incident_face, true)) {
    cout << "2" << endl;
    return true;
  }

  if (fix_proj(working_edge, projected, incident_face, prev_point,
               next_point)) {
    cout << "3" << endl;
    return true;
  }

  // tries to add triangle with prev point, true for prev
  if (fix_prev_next(working_edge, incident_face, true, true)) {
    cout << "4" << endl;
    return true;
  }
  // tries to add triangle with next point, false for next
  if (fix_prev_next(working_edge, incident_face, false, true)) {
    cout << "5" << endl;
    return true;
  }

  my_mesh.edges_check("inside step:");
  assertm(!is_border(working_edge_index),
          "Working edge found in border edges!");
  push_edge_to_checked(working_edge_index);
  my_mesh.edges_check("after step:");
  return false;
}
/*
bool BasicAlgorithm::fix_holes2(const Edge &working_edge) {
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);

  // if there is a hole in triangle shape, fill it with triangle
  if (basic_triangle(working_edge, neighbour_triangle, prev, next)) {
    return true;
  }
  // find candidate point for working_edge
  Point projected = get_projected(working_edge, neighbour_triangle);
  std::optional<Point> closest_point =
      get_closest_point(working_edge, neighbour_triangle);
  if (closest_point.has_value()) {
    create_triangle(working_edge, closest_point.value());
    return true;
  }

  if (fix_breakers(working_edge, projected, neighbour_triangle, false)) {
    return true;
  }
  if (fix_prev_next(working_edge, neighbour_triangle, true, false)) {
    return true;
  }
  if (fix_prev_next(working_edge, neighbour_triangle, false, false)) {
    return true;
  }

  push_edge_to_checked(working_edge);
  return false;
}
*/
void BasicAlgorithm::starting() {
  int round = 0;
  // cout << "in starting" << endl;
  while (!my_mesh.active_edges_empty()) {
    // cout << "round " << round << endl;
    round++;
    HalfEdgeIndex working_edge = kInvalidEdgeIndex;

    working_edge = my_mesh.get_active_edge();
    // my_mesh.obj_format(name);
    // my_mesh.cout_triangles();
    assertm(!my_mesh.is_in_mesh(
                Edge(my_mesh.get_halfedge(working_edge).get_point_B(),
                     my_mesh.get_halfedge(working_edge).get_point_A())),
            "Opposite edge found in mesh esges!");
    my_mesh.remove_active_edge(working_edge);

    assertm(working_edge != kInvalidEdgeIndex, "Active edge is invalid!");
    // cout << "before step" << endl;
    //  step returns true if new triangle is created
    if (step(working_edge))
      assertm(!my_mesh.is_boundary(working_edge),
              "Working edge found in border!");
    else
      assertm(my_mesh.is_checked(working_edge), "Checked edge not checked");

    // once in every 10 steps prints number of triangles and number of active
    // edges and updates output file
    // cout << "after step" << endl;
    my_mesh.obj_format(name);
    // my_mesh.cout_triangles_number();
    if (round % 50 == 0) {
      cout << "Number of active edges: " << my_mesh.get_active_edges_size()
           << endl;
      cout << "Number of checked edges: " << my_mesh.get_checked_edges_size()
           << endl
           << endl;
    }
  }

  // output
  my_mesh.obj_format(name);

  /*
  // this means there are no holes so we finished
  if (checked_edges.empty()) {
    cout << "No holes!" << endl;
    return;
  }
  return;
  */
}
/*
void BasicAlgorithm::ending() {
  assertm(active_edges.empty(), "Called ending with non-empty active edges!");

  // else push checked edges to active, clear checked and call ending

  int round0 = 0;
  while (!checked_edges.empty() && round0 < 2) {
    round0++;
    active_edges = checked_edges;
    checked_edges.clear();

    int round = 0;
    cout << "Fixing holes!" << endl;
    while (!active_edges.empty()) {
      round++;
      std::optional<Edge> working_edge = std::nullopt;

      working_edge = active_edges.back();
      active_edges.pop_back();

      assertm(working_edge.has_value(), "No working edge!");

      Triangle neighbour_triangle =
          my_mesh.find_triangle_with_edge(working_edge.value());
      fix_holes2(working_edge.value());
      if (round % 50 == 0) {
        my_mesh.cout_triangles_number();
        cout << "Number of active edges: " << active_edges.size() << endl
             << endl;
      }
    }
  }
  if (!checked_edges.empty()) {
    cout << "Some holes stayed!" << endl;
    for (auto edge : checked_edges) {
      cout << edge << endl;
    }
  }
  return;
}
*/

Mesh BasicAlgorithm::calculate() {
  my_mesh.edges_check("before starting:");
  starting();
  // ending();
  // fix_corners();
  //  my_mesh.measure(bounding_edges, F, name, e_size);
  //  my_mesh.adaptive(0.005, F, e_size);

  my_mesh.obj_format(name);
  return std::move(my_mesh);
}

// private member functions

// finds neighbour of prev which has the smallest angle with the working edge
std::optional<HalfEdgeIndex> BasicAlgorithm::_find_closest_prev(
    const HalfEdge &working_edge, const Face &incident_face,
    const vector<HalfEdgeIndex> &prev) const {
  std::optional<numeric> min_prev_angle = std::nullopt;
  std::optional<HalfEdgeIndex> min_prev_halfedge = std::nullopt;

  std::optional<numeric> my_angle = std::nullopt;

  for (auto prev_halfedge : prev) {
    assertm(prev_halfedge != working_edge.get_index(),
            "Prev halfedge is working edge!");

    // we will use angle function to find smallest angle
    my_angle = angle(my_mesh, working_edge,
                     my_mesh.get_halfedge(prev_halfedge).get_point_A(),
                     incident_face, false);
    if (!my_angle.has_value()) continue;
    if (my_angle.value() < 0)
      my_angle = my_angle.value() + ex_to<numeric>(2 * Pi.evalf());
    assertm(my_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                my_angle.value() >= 0,
            "Wrong angle interval!");

    if (min_prev_angle.has_value()) {
      if (my_angle.value() < min_prev_angle.value()) {
        min_prev_angle = my_angle.value();
        min_prev_halfedge = prev_halfedge;
      }
    } else {
      min_prev_angle = my_angle.value();
      min_prev_halfedge = prev_halfedge;
    }
    assertm(min_prev_angle.has_value(), "Angle without value!");
    assertm(min_prev_halfedge.has_value(), "Point without value!");
    assertm(min_prev_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                min_prev_angle.value() >= 0,
            "Angle not in right interval!");
  }

  return min_prev_halfedge;
}

// finds neighbour of next which has the smallest angle with the working edge
std::optional<HalfEdgeIndex> BasicAlgorithm::_find_closest_next(
    const HalfEdge &working_edge, const Face &incident_face,
    const vector<HalfEdgeIndex> &next) const {
  std::optional<numeric> min_next_angle = std::nullopt;
  std::optional<HalfEdgeIndex> min_next_halfedge = std::nullopt;

  std::optional<numeric> my_angle = std::nullopt;

  for (HalfEdgeIndex next_halfedge : next) {
    assertm(next_halfedge != working_edge.get_index(),
            "Working edge found amongst next halfedges!");
    const Point next_point = my_mesh.get_halfedge(next_halfedge).get_point_B();
    assertm(next_point != working_edge.get_point_B(),
            "Wrong potential next point B!");
    assertm(!my_mesh.is_in_mesh(
                Edge(working_edge.get_point_B(), working_edge.get_point_A())),
            "Opposite edge found in mesh esges!");
    assertm(next_point != working_edge.get_point_A(),
            "Wrong potential next point A!");
    // TODO invalid triangle still can have angle - 180
    /*assertm(Triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                     next_point)
                .is_triangle(),
            "Invalid triangle!");*/

    my_angle = angle(my_mesh, working_edge,
                     my_mesh.get_halfedge(next_halfedge).get_point_B(),
                     incident_face, true);
    if (!my_angle.has_value()) continue;

    if (my_angle.value() < 0)
      my_angle = my_angle.value() + ex_to<numeric>(2 * Pi.evalf());
    assertm(my_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                my_angle.value() >= 0,
            "Wrong angle interval!");

    if (min_next_angle.has_value()) {
      if (my_angle.value() < min_next_angle.value()) {
        min_next_angle = my_angle.value();
        min_next_halfedge = next_halfedge;
      }
    } else {
      min_next_angle = my_angle.value();
      min_next_halfedge = next_halfedge;
    }

    assertm(min_next_angle.has_value(), "Angle without value!");
    assertm(min_next_halfedge.has_value(), "Point without value!");
    assertm(min_next_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                min_next_angle.value() >= 0,
            "Angle not in right interval!");
  }

  return min_next_halfedge;
}

HalfEdgeIndex BasicAlgorithm::_find_next(const HalfEdge &working_edge,
                                         const Face &incident_face) const {
  const MeshPoint &meshpoint_next = my_mesh.get_meshpoint(working_edge.get_B());
  vector<HalfEdgeIndex> next;

  for (HalfEdgeIndex next_index : meshpoint_next.get_outgoing()) {
    if (my_mesh.get_halfedge(next_index).is_boundary()) {
      next.push_back(next_index);
    }
  }
  assertm(!next.empty(), "Neighbour edge not found in border edges!");

  if (next.size() == 1) return next[0];

  auto closest_next = _find_closest_next(working_edge, incident_face, next);

  assertm(closest_next.has_value(), "Closest next vithout value!");

  return closest_next.value();
}

HalfEdgeIndex BasicAlgorithm::_find_prev(const HalfEdge &working_edge,
                                         const Face &incident_face) const {
  const MeshPoint &meshpoint_prev = my_mesh.get_meshpoint(working_edge.get_A());

  vector<HalfEdgeIndex> prev_opposite = meshpoint_prev.get_outgoing();
  vector<HalfEdgeIndex> prev;

  // previous of edge outgoing from meshpoint_prev is incoming to meshpoint_prev
  for (HalfEdgeIndex prev_opp : prev_opposite) {
    // it should be previous, dont panic!
    HalfEdgeIndex prev_edge = my_mesh.get_previous_index(prev_opp);
    if (my_mesh.is_boundary(prev_edge)) prev.push_back(prev_edge);
  }

  assertm(!prev.empty(), "Neighbour edge not found in border edges!");

  if (prev.size() == 1) return prev[0];

  std::optional<HalfEdgeIndex> closest_prev =
      _find_closest_prev(working_edge, incident_face, prev);

  assertm(closest_prev.has_value(), "Closest prev or next without value!");

  return closest_prev.value();
}
