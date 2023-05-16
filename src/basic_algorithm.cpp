#include "basic_algorithm.h"

#include <iostream>

using std::cout;

#pragma region "Region singularities"

// given the polyline approximation of the singular curve
// creates the local mesh around the curve
// type - intersection, union, difference
void BasicAlgorithm::get_local_mesh_curve(
    const Function &_F, const Function &_G, const Function &inter_FG,
    const vector<Point> &polyline, const BoundingBox &bounding_box,
    const numeric e_size, const string type) {
  assertm(polyline.size() > 1, "Wrong length of polyline!");

  const Function *F = &_F;
  const Function *G = &_G;

  int offset_index_points = my_mesh.get_points_count();
  int offset_index_edges = my_mesh.get_edges_count();

  for (int i = 1; i < polyline.size(); ++i) {
    const Point &A = polyline[i - 1];
    const Point &B = polyline[i];
    const Point mid = (1.0 / 2.0) * (A + B);
    const std::optional<Point> opt_proj_F =
        project(mid, F->get_gradient_at_point(mid).unit(), *F, e_size);
    assertm(opt_proj_F.has_value(), "No value!");
    const Point proj_F = opt_proj_F.value();
    const std::optional<Point> opt_proj_G =
        project(mid, G->get_gradient_at_point(mid).unit(), *G, e_size);
    assertm(opt_proj_G.has_value(), "No value!");
    const Point proj_G = opt_proj_G.value();
    const Vector grad_F = F->get_gradient_at_point(proj_F).unit();
    const Vector grad_G = G->get_gradient_at_point(proj_G).unit();
    const Vector AB(B.x() - A.x(), B.y() - A.y(), B.z() - A.z());
    Vector tangent_F = e_size * (grad_F ^ AB).unit();
    Vector tangent_G = e_size * (grad_G ^ AB).unit();
    const Point tan_F_point(mid.x() + 0.2 * tangent_F.x(),
                            mid.y() + 0.2 * tangent_F.y(),
                            mid.z() + 0.2 * tangent_F.z());
    const Point tan_G_point(mid.x() + 0.2 * tangent_G.x(),
                            mid.y() + 0.2 * tangent_G.y(),
                            mid.z() + 0.2 * tangent_G.z());
    if (type == "intersection") {
      if (F->is_outside(tan_G_point)) tangent_G = -1 * tangent_G;
      if (G->is_outside(tan_F_point)) tangent_F = -1 * tangent_F;
    } else if (type == "union") {
      if (F->is_inside(tan_G_point)) tangent_G = -1 * tangent_G;
      if (G->is_inside(tan_F_point)) tangent_F = -1 * tangent_F;
    } else if (type == "difference") {
      if (F->is_outside(tan_G_point)) tangent_G = -1 * tangent_G;
      if (G->is_inside(tan_F_point)) tangent_F = -1 * tangent_F;
    } else {
      assertm(false, "Invalid type!");
    }
    const Point F_to_project(mid.x() + tangent_F.x(), mid.y() + tangent_F.y(),
                             mid.z() + tangent_F.z());
    const Point G_to_project(mid.x() + tangent_G.x(), mid.y() + tangent_G.y(),
                             mid.z() + tangent_G.z());
    std::optional<Point> opt_projected_F =
        project(F_to_project, F->get_gradient_at_point(F_to_project).unit(),
                inter_FG, e_size);
    assertm(opt_projected_F.has_value(), "No value!");
    Point projected_F = opt_projected_F.value();
    projected_F =
        bounding_box.crop_to_box(proj_F, projected_F, e_size, inter_FG);
    std::optional<Point> opt_projected_G =
        project(G_to_project, G->get_gradient_at_point(G_to_project).unit(),
                inter_FG, e_size);
    assertm(opt_projected_G.has_value(), "No value!");
    Point projected_G = opt_projected_G.value();
    projected_G =
        bounding_box.crop_to_box(proj_G, projected_G, e_size, inter_FG);
    const Triangle test_T(A, B, projected_F);

    /*
        if (test_T.get_normal() * inter_FG.outside_normal(test_T, e_size) < 0) {
          std::swap(projected_F, projected_G);
          std::swap(F, G);
        }
        */

    const Triangle T_F(A, B, projected_F);
    const Triangle T_G(B, A, projected_G);
    if (i == 1) {
      my_mesh.add_first_triangle(T_F, bounding_box);
    } else if (i == polyline.size() - 1 && polyline[0] == polyline[i]) {
      // special case - closed curve
      my_mesh.add_triangle_to_two_meshpoints(
          offset_index_points + 1 + (i - 2) * 3, offset_index_points + 0,
          projected_F, bounding_box);
    } else {
      my_mesh.add_triangle_to_meshpoint(offset_index_points + 1 + (i - 2) * 3,
                                        B, projected_F, bounding_box);
    }

    my_mesh.add_triangle(offset_index_edges + (i - 1) * 6, projected_G, true,
                         bounding_box);
  }
  // my_mesh.obj_format(name + "_local");
  return;
}

// rotation of point about vector by given angle
void rotate(Vector &to_rotate, const Vector &u, const GiNaC::ex angle) {
  const Vector v = u.get_any_perpendicular().unit();
  const numeric COS = GiNaC::ex_to<numeric>(GiNaC::cos(angle).evalf());
  const numeric SIN = GiNaC::ex_to<numeric>(GiNaC::sin(angle).evalf());
  const Vector rot_vector_x = Vector(COS + u.x() * u.x() * (1 - COS),
                                     u.x() * u.y() * (1 - COS) - u.z() * SIN,
                                     u.x() * u.z() * (1 - COS) + u.y() * SIN);
  const Vector rot_vector_y = Vector(u.x() * u.y() * (1 - COS) + u.z() * SIN,
                                     COS + u.y() * u.y() * (1 - COS),
                                     u.y() * u.z() * (1 - COS) - u.x() * SIN);
  const Vector rot_vector_z = Vector(u.x() * u.z() * (1 - COS) - u.y() * SIN,
                                     u.y() * u.z() * (1 - COS) + u.x() * SIN,
                                     COS + u.z() * u.z() * (1 - COS));
  to_rotate = Vector(rot_vector_x * to_rotate, rot_vector_y * to_rotate,
                     rot_vector_z * to_rotate);
  return;
}

// creates layers of points around the given branch of a singularity
void BasicAlgorithm::get_points_singular_point(const Singularity singularity,
                                               const int num_triangles,
                                               const numeric &edge_length,
                                               vector<vector<Point>> &points,
                                               const int branch,
                                               const int layers) {
  const Vector unit_singular_direction =
      singularity.get_direction(branch).unit();
  Vector plane_vector = unit_singular_direction.get_any_perpendicular().unit();
  const Vector &u = unit_singular_direction;
  const numeric COS = GiNaC::ex_to<numeric>(
      GiNaC::cos(2 * GiNaC::Pi.evalf() / num_triangles).evalf());
  const numeric SIN = GiNaC::ex_to<numeric>(
      GiNaC::sin(2 * GiNaC::Pi.evalf() / num_triangles).evalf());
  const Vector rot_vector_x = Vector(COS + u.x() * u.x() * (1 - COS),
                                     u.x() * u.y() * (1 - COS) - u.z() * SIN,
                                     u.x() * u.z() * (1 - COS) + u.y() * SIN);
  const Vector rot_vector_y = Vector(u.x() * u.y() * (1 - COS) + u.z() * SIN,
                                     COS + u.y() * u.y() * (1 - COS),
                                     u.y() * u.z() * (1 - COS) - u.x() * SIN);
  const Vector rot_vector_z = Vector(u.x() * u.z() * (1 - COS) - u.y() * SIN,
                                     u.y() * u.z() * (1 - COS) + u.x() * SIN,
                                     COS + u.z() * u.z() * (1 - COS));
  const numeric e = edge_length;
  numeric h;
  if (singularity.type() == SingularityType::Amm) {
    h = pow(e, (numeric)2 / (singularity.n() + 1));
  } else {
    h = edge_length;
  }
  numeric h_layer, e_layer;
  const numeric delta_angle = GiNaC::ex_to<numeric>(GiNaC::Pi.evalf() / 36);
  for (int layer = 0; layer < layers; ++layer) {
    plane_vector = unit_singular_direction.get_any_perpendicular().unit();
    h_layer = (layer + 1) * h / (numeric)layers;
    e_layer = pow(h_layer, ((numeric)singularity.n() + 1) / (numeric)2);
    Vector direction_layer = plane_vector * e_layer + u * h_layer;
    if (layer == 1) {
      rotate(direction_layer, u, GiNaC::Pi.evalf() / num_triangles);
      // rotate(plane_vector, u, GiNaC::Pi.evalf() / num_triangles);
    } else if (layer > 2) {
      rotate(direction_layer, u,
             -(layer - 2) * GiNaC::Pi.evalf() / num_triangles);
      // rotate(plane_vector, u, -(layer - 2) * GiNaC::Pi.evalf() /
      // num_triangles);
    }
    points.push_back(vector<Point>());
    for (int i = 0; i < num_triangles; ++i) {
      Point new_point;
      if (singularity.type() == SingularityType::Amm) {  // group 0 - analytical
        new_point = Point(singularity.x() + direction_layer.x(),
                          singularity.y() + direction_layer.y(),
                          singularity.z() + direction_layer.z());
      } else if ((singularity.type() == SingularityType::Apm) ||
                 (singularity.type() == SingularityType::Dpm &&
                  singularity.n() % 2 == 0) ||
                 (singularity.type() == SingularityType::Epp &&
                  singularity.n() % 2 == 0) ||
                 singularity.type() ==
                     SingularityType::Epm) {  // group 1 - circular bisection
        Point start_point_to_bisect, end_point_to_bisect;
        numeric angle;
        if (singularity.type() == SingularityType::Apm &&
            singularity.n() % 2 == 1) {
          start_point_to_bisect =
              Point(singularity.location(), h_layer * plane_vector);
          direction_layer = (plane_vector ^ unit_singular_direction).unit();
          angle = GiNaC::ex_to<numeric>((GiNaC::Pi / 2).evalf());
        } else {
          start_point_to_bisect =
              Point(singularity.location(), -h_layer * unit_singular_direction);

          direction_layer = (plane_vector).unit();
          angle = GiNaC::ex_to<numeric>((GiNaC::Pi).evalf());
        }
        end_point_to_bisect =
            Point(singularity.location(), h_layer * unit_singular_direction);
        new_point =
            circular_bisection(F, singularity.location(), start_point_to_bisect,
                               angle, direction_layer, h_layer);

        plane_vector =
            Vector(rot_vector_x * plane_vector, rot_vector_y * plane_vector,
                   rot_vector_z * plane_vector);
      } else {  // group 2 - iterative binary search
        Point start_point_to_bisect(singularity.location(),
                                    h_layer * unit_singular_direction);
        Point end_point_to_bisect = start_point_to_bisect;
        int iter = 0;
        while (F.eval_at_point(start_point_to_bisect) *
                   F.eval_at_point(end_point_to_bisect) >
               0) {
          assertm(iter <= 71, "Too many iterations!");
          end_point_to_bisect =
              rotate(singularity.location(), end_point_to_bisect, plane_vector,
                     delta_angle);
          iter++;
        }
        new_point =
            circular_bisection(F, singularity.location(), start_point_to_bisect,
                               iter * delta_angle, plane_vector, h_layer);
        plane_vector = Vector(
            singularity.location(),
            rotate(
                singularity.location(),
                Point(plane_vector.x() - singularity.x(),
                      plane_vector.y() - singularity.y(),
                      plane_vector.z() - singularity.z()),
                u,
                GiNaC::ex_to<numeric>(2 * GiNaC::Pi.evalf() / num_triangles)));
      }
      points[layer].push_back(new_point);

      // checks if the point is in the correct halfspace
      assertm(unit_singular_direction * direction_layer >= 0,
              "New point in opposite halfspace!");
      direction_layer =
          Vector(rot_vector_x * direction_layer, rot_vector_y * direction_layer,
                 rot_vector_z * direction_layer);
    }
    
    // the optimization of the mesh for some singularities with 8 triangles
    if (num_triangles == 8) {
      if (singularity.type() == SingularityType::Epp &&
          singularity.n() == 8)  // E8++ not symmetrical
        return;
      if (singularity.type() == SingularityType::Apm &&
          singularity.n() % 2 == 0)  // An+- not necessary
        return;
      numeric dist01 = Vector(points[layer][0], points[layer][1]).get_length();
      numeric dist12 = Vector(points[layer][1], points[layer][2]).get_length();
      numeric dist02 = Vector(points[layer][0], points[layer][2]).get_length();
      numeric dist_difference = abs(dist01 - dist12);
      numeric angle1 = 0;
      numeric angle2 = GiNaC::ex_to<numeric>((GiNaC::Pi / 2).evalf());
      plane_vector = unit_singular_direction.get_any_perpendicular().unit();
      Point new_point;
      Vector rotated_plane_vector = plane_vector;
      int iter = 0;
      while (dist_difference > dist02 / 20 && iter < 100) {
        iter++;
        numeric angle_mid = (angle1 + angle2) / 2;

        if ((singularity.type() == SingularityType::Dpm &&
             singularity.n() % 2 == 1) ||
            (singularity.type() == SingularityType::Dmm &&
             singularity.n() % 2 == 0)) {
          rotated_plane_vector =
              Vector(singularity.location(),
                     rotate(singularity.location(),
                            Point(plane_vector.x() - singularity.x(),
                                  plane_vector.y() - singularity.y(),
                                  plane_vector.z() - singularity.z()),
                            u, angle_mid));
          Point start_point_to_bisect(singularity.location(),
                                      h_layer * unit_singular_direction);
          Point end_point_to_bisect = start_point_to_bisect;
          int iter = 0;
          while (F.eval_at_point(start_point_to_bisect) *
                     F.eval_at_point(end_point_to_bisect) >
                 0) {
            assertm(iter <= 71, "Too many iterations!");
            end_point_to_bisect =
                rotate(singularity.location(), end_point_to_bisect,
                       rotated_plane_vector, delta_angle);
            iter++;
          }
          new_point = circular_bisection(
              F, singularity.location(), start_point_to_bisect,
              iter * delta_angle, rotated_plane_vector, h_layer);
        } else if ((singularity.type() == SingularityType::Apm) ||
                   (singularity.type() == SingularityType::Dpm &&
                    singularity.n() % 2 == 0) ||
                   (singularity.type() == SingularityType::Epp &&
                    singularity.n() % 2 == 0) ||
                   singularity.type() == SingularityType::Epm) {
          rotated_plane_vector =
              Vector(singularity.location(),
                     rotate(singularity.location(),
                            Point(plane_vector.x() - singularity.x(),
                                  plane_vector.y() - singularity.y(),
                                  plane_vector.z() - singularity.z()),
                            u, angle_mid));

          Point start_point_to_bisect, end_point_to_bisect;
          numeric angle;
          if (singularity.type() == SingularityType::Apm &&
              singularity.n() % 2 == 1) {
            start_point_to_bisect =
                Point(singularity.location(), h_layer * rotated_plane_vector);
            direction_layer =
                (rotated_plane_vector ^ unit_singular_direction).unit();
            angle = GiNaC::ex_to<numeric>((GiNaC::Pi / 2).evalf());
          } else {
            start_point_to_bisect = Point(singularity.location(),
                                          -h_layer * unit_singular_direction);

            direction_layer = (rotated_plane_vector).unit();
            angle = GiNaC::ex_to<numeric>((GiNaC::Pi).evalf());
          }
          end_point_to_bisect =
              Point(singularity.location(), h_layer * unit_singular_direction);
          new_point = circular_bisection(F, singularity.location(),
                                         start_point_to_bisect, angle,
                                         direction_layer, h_layer);
        }
        points[layer][1] = new_point;
        dist01 = Vector(points[layer][0], points[layer][1]).get_length();
        dist12 = Vector(points[layer][1], points[layer][2]).get_length();
        dist_difference = abs(dist01 - dist12);
        if (dist01 > dist12) {
          angle2 = angle_mid;
        } else {  // if(dist01 < dist12)
          angle1 = angle_mid;
        }
      }
      int first, second;
      const Point &P = points[layer][1];
      
      // TODO fix for general case, not working properly
      if (abs(u.x()) < kEps && abs(u.z()) < kEps) {
        first = 3;
        second = 7;
        points[layer][first] = Point(P.x(), P.y(), -P.z());
        points[layer][5] = Point(-P.x(), P.y(), -P.z());
        points[layer][second] = Point(-P.x(), P.y(), P.z());
      } else if (abs(u.y()) < kEps && abs(u.z()) < kEps) {
        first = 3;
        second = 7;
        points[layer][first] = Point(P.x(), P.y(), -P.z());
        points[layer][5] = Point(P.x(), -P.y(), -P.z());
        points[layer][second] = Point(P.x(), -P.y(), P.z());
      } else if (abs(u.x()) < kEps && abs(u.y()) < kEps) {
        first = 3;
        second = 7;
        points[layer][first] = Point(P.x(), -P.y(), P.z());
        points[layer][5] = Point(-P.x(), -P.y(), P.z());
        points[layer][second] = Point(-P.x(), P.y(), P.z());
      }
    }
  }
  return;
}

// connects the layers of points into triangles
void BasicAlgorithm::get_local_mesh_point(const vector<vector<Point>> &points,
                                          const Singularity &singularity,
                                          const MeshPointIndex singular_index,
                                          const int num_triangles,
                                          const int branch, const int layers) {
  // the patterns of in indices were derived on a paper

  // first layer
  for (int i = 0; i < num_triangles; ++i) {
    int j = (i + 1) % num_triangles;
    Triangle triangle =
        Triangle(singularity.location(), points[0][i], points[0][j]);
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (my_mesh.get_faces_count() == 0 ||
          singular_index == kInvalidPointIndex) {
        my_mesh.add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        my_mesh.add_triangle_to_meshpoint(singular_index, points[0][i],
                                          points[0][j], bounding_box);
      }
    } else if (i < num_triangles - 1) {
      my_mesh.add_triangle(my_mesh.get_edges_count() - 1, points[0][j], true,
                           bounding_box);
    } else {
      my_mesh.add_triangle(my_mesh.get_edges_count() - 1, points[0][j], false,
                           bounding_box,
                           my_mesh.get_points_count() - num_triangles);
    }
  }

  // the rest of the layers
  for (int layer = 1; layer < layers; ++layer) {
    int magic_constant = layer == 1 ? 1 : 2;
    for (int i = 0; i < num_triangles; ++i) {
      my_mesh.add_triangle(
          my_mesh.get_faces_count() * 3 - (3 * num_triangles - magic_constant),
          points[layer][i], true, bounding_box);
    }
    for (int i = 0; i < num_triangles; ++i) {
      int k = (i + num_triangles - 1) % num_triangles;
      my_mesh.add_triangle(
          3 * my_mesh.get_faces_count() - (3 * num_triangles - 1),
          points[layer][k], false, bounding_box,
          my_mesh.get_points_count() - num_triangles + k);
    }
  }
}

// division into cases - possible twisting of parameters to some extent
void BasicAlgorithm::triangulate_singular_point_local(
    const Singularity &singularity, const int branch,
    const MeshPointIndex singular_index) {
  int num_triangles;
  int layers;
  numeric length;
  if (singularity.type() == SingularityType::Amm) {  // An--
    num_triangles = 6;
    layers = singularity.n() + 1;
    length = e_size;
  } else if (singularity.type() == SingularityType::Apm &&
             singularity.n() % 2 == 0)  // A2+-, A4+-, A6+-, ...
  {
    num_triangles = 8;
    layers = 1;  // singularity.n()+1;
    length = e_size;
  } else if (singularity.type() == SingularityType::Apm &&
             singularity.n() % 2 == 1) {  // A1+-, A3+-, A5+-, ...
    num_triangles = 6;
    layers = 1;
    length = e_size;
  } else if (singularity.type() == SingularityType::Dpm &&
             singularity.n() % 2 == 0) {  // D4+-, D6+-, D8+-, ...
    num_triangles = 8;
    layers = 2;
    length = e_size;
  } else if (singularity.type() == SingularityType::Epm) {  // E6+-
    num_triangles = 8;
    layers = 2;
    length = e_size * 1.5;
  } else if (singularity.type() == SingularityType::Epp &&
             singularity.n() % 2 == 0) {  // E6++, E8++
    num_triangles = 4;                    // 6, 8;
    layers = 1;
    length = e_size * 1.5;
  } else if (singularity.type() == SingularityType::Epp &&
             singularity.n() % 2 == 1) {  // E7++
    if (branch == 0) {
      num_triangles = 6;
      layers = 2;
      length = e_size;
    } else {
      num_triangles = 7;  // 4, 8;
      layers = 3;
      length = layers * e_size;
    }
  } else if (singularity.type() == SingularityType::Dpm &&
             singularity.n() % 2 == 1) {  // D5+-, D7+-, D9+-, ...
    if (branch == 0) {
      num_triangles = 8;
      layers = 4;
      length = 2 * e_size;
    } else {
      num_triangles = 8;
      layers = 4;
      length = 2 * e_size;
    }
  } else if (singularity.type() == SingularityType::Dmm &&
             singularity.n() % 2 == 0) {  // D4--, D6--, ...
    if (branch == 0) {
      num_triangles = 8;
    } else {
      num_triangles = 7;
    }
    layers = 3;
    length = layers * e_size;
  } else if (singularity.type() == SingularityType::Dmm &&
             singularity.n() % 2 == 1) {  // D5--, D7--, ...
    num_triangles = 7;

    layers = 3;
    length = layers * e_size;
  } else {
    assertm(false, "Should not get here!");
  }
  vector<vector<Point>> points;
  get_points_singular_point(singularity, num_triangles, length, points, branch,
                            layers);
  get_local_mesh_point(points, singularity, singular_index, num_triangles,
                       branch, layers);
  return;
}
#pragma endregion "Region singularities"

#pragma region "Find close points"

// using the range tree to find the points in the surrounding of a point P
// sorted by the distance from the edge introduced in bachelor's thesis
std::optional<vector<MeshPoint>> BasicAlgorithm::_tree_close_points_finder(
    const Point &P, const HalfEdge &working_edge) const {
  numeric dist = 0.4 * e_size;
  std::vector<MeshPoint> close_points = my_mesh.get_meshpoints_in_interval(
      P.x() - dist, P.x() + dist, P.y() - dist, P.y() + dist, P.z() - dist,
      P.z() + dist);

  if (close_points.empty()) return std::nullopt;
  sort(close_points.begin(), close_points.end(),
       [this, &working_edge](const auto &i, const auto &j) {
         return my_mesh.line_point_dist(working_edge, i) <
                my_mesh.line_point_dist(working_edge, j);
       });
  return close_points;
}

// using the linear search to find the points in the surrounding of a point P
// sorted by the distance from the edge introduced in bachelor's thesis
std::optional<vector<MeshPoint>> BasicAlgorithm::_linear_close_points_finder(
    const Point &P, const HalfEdge &working_edge) const {
  numeric min_dist = 0.4 * e_size;
  vector<MeshPoint> close_points;
  for (MeshPoint meshpoint : my_mesh.get_mesh_points()) {
    if (my_mesh.is_boundary_point(meshpoint) &&
        Vector(meshpoint.get_point(), P).get_length() < min_dist)
      close_points.push_back(meshpoint);
  }
  if (close_points.empty()) return std::nullopt;

  sort(close_points.begin(), close_points.end(),
       [this, &working_edge](const auto &i, const auto &j) {
         return my_mesh.line_point_dist(working_edge, i) <
                my_mesh.line_point_dist(working_edge, j);
       });
  return close_points;
}

// returns vector of points closer than 0.4*e_size to point P sorted from
// closest to working edge if there are no points returns std::nullopt
std::optional<vector<MeshPoint>> BasicAlgorithm::find_close_points(
    const Point &P, const HalfEdge &working_edge) const {
  return _tree_close_points_finder(P, working_edge);
}

#pragma endregion "Find close points"

#pragma region "Check conditions"

// check of the conditions for the creation of the new triangle
bool BasicAlgorithm::check_conditions(const HalfEdge &working_edge,
                                      const Point &P, const bool Delaunay,
                                      const MeshPointIndex i_P) const {
  if (Delaunay) return Delaunay_conditions(working_edge, P);
  return non_Delaunay_conditions(working_edge, P, i_P);
}

// checks if conditions required in the first part of the algorithm are
// satisfied
bool BasicAlgorithm::Delaunay_conditions(const HalfEdge &working_edge,
                                         const Point &P) const {
  if (working_edge.get_point_A() == P || working_edge.get_point_B() == P) {
    return false;
  }
  Triangle new_triangle =
      Triangle(working_edge.get_point_B(), working_edge.get_point_A(), P);

  if (!new_triangle.is_triangle()) {
    return false;
  }

  numeric min_edge_length = std::min(
      std::min(new_triangle.AB().get_length(), new_triangle.BC().get_length()),
      new_triangle.CA().get_length());
  numeric avg_edge_length =
      (new_triangle.AB().get_length() + new_triangle.BC().get_length() +
       new_triangle.CA().get_length()) /
      3.0;
  
  // possible twisting the parameter 2.5
  if (/*my_mesh.get_faces_count() > 50 &&*/ // turned on for some singularities with sharp triangles 
      (new_triangle.AB().get_length() > 2.5 * min_edge_length ||
       new_triangle.CA().get_length() > 2.5 * min_edge_length ||
       new_triangle.BC().get_length() > 2.5 * min_edge_length)) {
    return false;
  }

  const Point opposite_point =
      my_mesh
          .get_meshpoint(my_mesh.get_halfedge(working_edge.get_next()).get_B())
          .get_point();

  if (P == opposite_point) return false;

  // possible to turn on in case of the connection of two distinct branches
  bool gc_distance = true;  // gc_distance_check(new_triangle);
  bool delaunay = my_mesh.check_Delaunay(working_edge, P, avg_edge_length);
  bool has_good_edges = good_edges(working_edge, P);
  bool neighbor_triangles_normal_check = _check_normals(working_edge, P);
  bool is_edge_in_mesh = my_mesh.is_in_mesh(new_triangle.AB()) ||
                         my_mesh.is_in_mesh(new_triangle.BC()) ||
                         my_mesh.is_in_mesh(new_triangle.CA());
  
  // bool orientability = orientability_check(working_edge, P);
  return (delaunay && has_good_edges && neighbor_triangles_normal_check &&
          !is_edge_in_mesh && gc_distance
          // && orientability
  );
}

// checks if essential conditions in the hole fixing are fulfilled
bool BasicAlgorithm::non_Delaunay_conditions(const HalfEdge &working_edge,
                                             const Point &P,
                                             const MeshPointIndex i_P) const {
  if (working_edge.get_point_A() == P || working_edge.get_point_B() == P) {
    return false;
  }
  Triangle new_triangle =
      Triangle(working_edge.get_point_B(), working_edge.get_point_A(), P);

  if (!new_triangle.is_triangle()) return false;
  bool has_good_edges = good_edges(working_edge, P);
  assertm(i_P != kInvalidPointIndex, "Invalid meshpoint index!");
  bool is_border = my_mesh.is_boundary_point(my_mesh.get_meshpoint(i_P));
  bool is_edge_in_mesh = my_mesh.is_in_mesh(new_triangle.AB()) ||
                         my_mesh.is_in_mesh(new_triangle.BC()) ||
                         my_mesh.is_in_mesh(new_triangle.CA());
  // bool orientability = orientability_check(working_edge, P);

  return (has_good_edges && is_border && !is_edge_in_mesh  //&& orientability
  );
}

// we can only create a new tringle with boundary meshpoint
bool BasicAlgorithm::good_new_point(const MeshPoint &point) const {
  return my_mesh.is_boundary_point(point);
}

// checks if edges of new triangle are active or are not im mesh
bool BasicAlgorithm::good_edges(const HalfEdge &working_edge,
                                const Point &P) const {
  if (working_edge.get_point_A() == P || working_edge.get_point_B() == P)
    return false;
  Edge new_edge1(P, working_edge.get_point_A());
  Edge new_edge2(working_edge.get_point_B(), P);

  HalfEdgeIndex new_edge1_index = my_mesh.get_edge_index(new_edge1);
  HalfEdgeIndex new_edge2_index = my_mesh.get_edge_index(new_edge2);

  // Edges are good if either they are not in mesh
  // or they are in mesh and are border

  // TODO(kuska) add good egdes condition and check perpendicular to edge
  // if new_edge1 is in mesh it can be only incoming to A
  // if new_edge2 is in mesh it can be only outgoing of B

  bool is_in_mesh_1 = (new_edge1_index != kInvalidEdgeIndex);
  bool is_in_mesh_2 = (new_edge2_index != kInvalidEdgeIndex);

  if (!is_in_mesh_1 && !is_in_mesh_2) return true;

  bool is_good_1 = true, is_good_2 = true;

  if (is_in_mesh_1) is_good_1 = is_border(new_edge1_index);
  if (is_in_mesh_2) is_good_2 = is_border(new_edge2_index);

  return (is_good_1 && is_good_2);
}

bool BasicAlgorithm::gc_distance_check(const Triangle T) const {
  const Point gravity_center = T.get_gravity_center();
  const numeric min_edge_length = std::min(
      T.AB().get_length(), std::min(T.BC().get_length(), T.CA().get_length()));
  const std::optional<Point> proj_gravity_center =
      project(gravity_center, F.get_gradient_at_point(gravity_center).unit(), F,
              e_size);
  if (!proj_gravity_center.has_value()) return false;
  if (Vector(gravity_center, proj_gravity_center.value()).get_length() >
      min_edge_length / 10)
    return false;

  const Point M_AB = T.AB().get_midpoint();
  const Point M_BC = T.BC().get_midpoint();
  const Point M_CA = T.CA().get_midpoint();

  const std::optional<Point> proj_M_AB =
      project(M_AB, F.get_gradient_at_point(M_AB).unit(), F, e_size);
  const std::optional<Point> proj_M_BC =
      project(M_BC, F.get_gradient_at_point(M_BC).unit(), F, e_size);
  const std::optional<Point> proj_M_CA =
      project(M_CA, F.get_gradient_at_point(M_CA).unit(), F, e_size);
  if (!proj_M_AB.has_value()) return false;
  if (!proj_M_BC.has_value()) return false;
  if (!proj_M_CA.has_value()) return false;
  if (Vector(M_AB, proj_M_AB.value()).get_length() > T.AB().get_length() / 4)
    return false;
  if (Vector(M_BC, proj_M_BC.value()).get_length() > T.BC().get_length() / 4)
    return false;
  if (Vector(M_CA, proj_M_CA.value()).get_length() > T.CA().get_length() / 4)
    return false;

  return true;
}

bool BasicAlgorithm::neighbor_triangles_normal_check(const Triangle &T1,
                                                     const Triangle &T2) const {
  assertm(T1.is_triangle() && T2.is_triangle(),
          "Normal check of invalid triangle");

  Vector normal_1 = T1.get_normal();
  Vector normal_2 = T2.get_normal();
  return (normal_1 * normal_2 > 0.25);
}

// returns false if the new edges already exist in mesh
bool BasicAlgorithm::orientability_check(const HalfEdge &working_edge,
                                         const MeshPoint &point) const {
  const MeshPoint &A = my_mesh.get_meshpoint(working_edge.get_A());
  const MeshPoint &B = my_mesh.get_meshpoint(working_edge.get_B());
  // AP
  for (HalfEdgeIndex outgoing : A.get_outgoing()) {
    if (my_mesh.get_halfedge(outgoing).get_B() == point.get_index())
      return false;
  }
  // BA
  for (HalfEdgeIndex outgoing : B.get_outgoing()) {
    if (my_mesh.get_halfedge(outgoing).get_B() == working_edge.get_A())
      return false;
  }
  // PB
  for (HalfEdgeIndex outgoing : point.get_outgoing()) {
    if (my_mesh.get_halfedge(outgoing).get_B() == working_edge.get_B())
      return false;
  }
  return true;
}

// checks if the incident triangles have angle smaller than pi/2
bool BasicAlgorithm::overlap_normals_check(const MeshPoint &candidate,
                                           const HalfEdge &working_edge) const {
  if (candidate.get_point() == working_edge.get_point_A() ||
      candidate.get_point() == working_edge.get_point_B()) {
    return false;
  }

  Triangle my_triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                       candidate.get_point());

  if (my_triangle.is_triangle()) {
    Vector my_normal = my_triangle.get_normal();
    Face overlap_face = my_mesh.get_face(
        my_mesh.get_halfedge(candidate.get_outgoing()[0]).get_incident());
    Vector overlap_normal = overlap_face.get_normal();
    // if normals have the same orientation
    if (overlap_normal * my_normal > 0) {
      return true;
    }
    return false;
  }
  return false;
}

#pragma endregion "Check conditions"

// creates new triangle and adds it to mesh
void BasicAlgorithm::create_triangle(const HalfEdge &working_edge,
                                     const Point &P, const bool is_new,
                                     const MeshPointIndex index_P) {
  assertm(working_edge.get_point_A() != P && working_edge.get_point_B() != P,
          "Wrong new triangle!");
  Triangle new_triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                        P);
  assertm(new_triangle.is_triangle(), "New triangle is not a triangle!");
  my_mesh.add_triangle(working_edge.get_index(), P, is_new, bounding_box,
                       index_P);
  return;
}

// if there is a triangle hole, returns index of the third vertex
MeshPointIndex BasicAlgorithm::find_triangle_hole(
    const HalfEdge &working_edge) {
  const MeshPoint &A = my_mesh.get_meshpoint(working_edge.get_A());
  const MeshPoint &B = my_mesh.get_meshpoint(working_edge.get_B());

  const MeshPointIndex &opposite_point =
      my_mesh.get_halfedge(working_edge.get_next()).get_B();

  for (HalfEdgeIndex outgoing_B : B.get_outgoing()) {
    // if (!my_mesh.is_boundary(outgoing_B)) continue;
    const HalfEdge &outgoing_edge_B = my_mesh.get_halfedge(outgoing_B);
    const MeshPointIndex &outgoing_point_B = outgoing_edge_B.get_B();
    if (outgoing_point_B == opposite_point) continue;
    for (HalfEdgeIndex outgoing_A : A.get_outgoing()) {
      const HalfEdge &incoming_A = my_mesh.get_previous_halfedge(outgoing_A);
      // if (!my_mesh.is_boundary(incoming_A.get_index())) continue;
      const MeshPointIndex &incoming_point_A = incoming_A.get_A();
      if (incoming_point_A == outgoing_point_B) {
        return incoming_point_A;
      }
    }
  }
  return kInvalidPointIndex;
}

// finds incident faces of a potential triangle
std::pair<FaceIndex, FaceIndex> BasicAlgorithm::find_neighbor_faces(
    const HalfEdge &working_edge, const Point &P) const {
  const MeshPoint &A = my_mesh.get_meshpoint(working_edge.get_A());
  const MeshPoint &B = my_mesh.get_meshpoint(working_edge.get_B());

  FaceIndex previous_face = kInvalidFaceIndex;
  FaceIndex next_face = kInvalidFaceIndex;

  for (HalfEdgeIndex outgoing : B.get_outgoing()) {
    const HalfEdge &outgoing_edge = my_mesh.get_halfedge(outgoing);
    const MeshPoint &outgoing_point =
        my_mesh.get_meshpoint(outgoing_edge.get_B());
    if (P == outgoing_point.get_point()) {
      next_face = outgoing_edge.get_incident();
      break;
    }
  }
  for (HalfEdgeIndex outgoing : A.get_outgoing()) {
    const HalfEdge &outgoing_edge = my_mesh.get_halfedge(outgoing);
    const HalfEdge &incoming_edge =
        my_mesh.get_halfedge(outgoing_edge.get_previous());
    const MeshPoint &incoming_point =
        my_mesh.get_meshpoint(incoming_edge.get_A());
    if (P == incoming_point.get_point()) {
      previous_face = incoming_edge.get_incident();
      break;
    }
  }
  return {previous_face, next_face};
}

// fills in triangle hole
bool BasicAlgorithm::basic_triangle(const HalfEdge &working_edge,
                                    const MeshPoint &point) {
  // determine other point on the neighbour triangle - opposite to working_edge
  MeshPoint opposite_point = my_mesh.get_meshpoint(
      my_mesh.get_halfedge(working_edge.get_next()).get_B());
  assertm(my_mesh.has_incoming_prev(working_edge, point.get_point()) &&
              my_mesh.has_outgoing_next(working_edge, point.get_point()),
          "Wrong call of basic triangle!");
  if (working_edge.get_point_A() == point.get_point() ||
      working_edge.get_point_B() == point.get_point())
    return false;
  const Triangle new_triangle =
      Triangle(working_edge.get_point_B(), working_edge.get_point_A(),
               point.get_point());
  assertm(new_triangle.is_triangle(), "Wrong call of basic triangle!");

  // std::cout << "checking conditions for: " << working_edge <<
  // point.get_point()
  //           << endl;
  if (!check_conditions(working_edge, point.get_point(), false,
                        point.get_index())) {
    return false;
  }
  /*std::cout << "neighbour triangle: "
            << my_mesh.get_face(working_edge.get_incident()).get_triangle()
            << endl;
  std::cout << "new triangle: " << new_triangle << endl;*/
  if (opposite_point == point) return false;
  create_triangle(working_edge, point.get_point(), false, point.get_index());
  return true;
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
/*
std::optional<MeshPoint> BasicAlgorithm::get_closest_point(
    const HalfEdge &working_edge) const {
  std::vector<MeshPoint> meshpoints = my_mesh.get_mesh_points();
  std::optional<MeshPoint> closest_point = std::nullopt;
  std::optional<numeric> min_dist;

  for (auto meshpoint : meshpoints) {
    auto dist = my_mesh.line_point_dist(working_edge, meshpoint.get_point());
    if (dist < min_dist &&
        meshpoint.get_point() != working_edge.get_point_A() &&
        meshpoint.get_point() != working_edge.get_point_B() &&
        check_conditions(working_edge, meshpoint.get_point(), false,
                         meshpoint.get_index()))
      min_dist = dist;
    closest_point = meshpoint;
  }
  if (closest_point.has_value() && min_dist.has_value() &&
      min_dist.value() < 1.5 * e_size) {
    return closest_point.value();
  }
  return std::nullopt;
}
*/

std::optional<Point> BasicAlgorithm::get_projected(
    const HalfEdge &working_edge) const {
  Point center = working_edge.get_midpoint();
  assertm(Vector(working_edge.get_point_A(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6 &&
              Vector(working_edge.get_point_B(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6,
          "Wrong get_midpoint function!");

  // height of equilateral triangle based on e_size

  /*
    numeric average_edge_length =
         (1 / numeric(3)) *
         (Edge(working_edge.get_point_A(), prev_point.get_point()).get_length()
         +
          Edge(working_edge.get_point_B(), next_point.get_point()).get_length()
          + working_edge.get_length());
  */
  // height of equilateral triangle based on working_edge size
  // numeric height = working_edge.get_length() * sqrt(numeric(3)) / 2;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;

  // height of equilateral triangle based on neighbour edges size
  // numeric height = average_edge_length * sqrt(numeric(3)) / 2;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;

  // height of equilateral triangle based on neighbour edges size with influence
  // of e_size
  /*
  numeric height = 0.75 * (working_edge.get_length() * sqrt(numeric(3)) / 2) +
                   0.25 * basic_height;
  if (height < e_size / 3) height = 0.25 * (e_size / 3) + 0.75 * height;
  if (height > 2 * e_size) height = 2 * e_size;
  */
  //  height = gaussian_height;
  //   non adaptive height
  //   numeric height = gaussian_height;

  // section -- change the incident face plane to tangent plane

  const Face &incident_face = my_mesh.get_face(working_edge.get_incident());
  const Vector gradient_at_edge_midpoint = F.get_gradient_at_point(center);
  numeric height;
  numeric basic_height = e_size * sqrt(numeric(3)) / 2;

  bool adaptive = false;
  if (adaptive && !gradient_at_edge_midpoint.is_zero()) {
    Vector normal_at_edge_midpoint = gradient_at_edge_midpoint.unit();
    if (normal_at_edge_midpoint * incident_face.get_normal() < 0)
      normal_at_edge_midpoint = -1 * normal_at_edge_midpoint;
    numeric edge_height = working_edge.get_length() * sqrt(numeric(3)) / 2;
    const Vector dir = find_direction_plane(
        working_edge, normal_at_edge_midpoint, incident_face);
    const Point midpoint = Point(center, (edge_height)*dir);

    numeric avg_e_size = working_edge.get_length();

    vector<MeshPoint> prev = my_mesh.get_prev(working_edge);
    vector<MeshPoint> next = my_mesh.get_next(working_edge);
    for (MeshPoint p : prev) {
      avg_e_size +=
          Vector(p.get_point(), working_edge.get_point_A()).get_length();
    }
    for (MeshPoint n : next) {
      avg_e_size +=
          Vector(n.get_point(), working_edge.get_point_B()).get_length();
    }
    avg_e_size /= (1 + prev.size() + next.size());

    std::optional<Point> proj_midpoint =
        project(center, F.get_gradient_at_point(midpoint).unit(), F, e_size);
    if (proj_midpoint.has_value()) {
      numeric gaussian_mult = 10000;
      if (!F.get_gradient_at_point(working_edge.get_point_A()).is_zero()) {
        gaussian_mult =
            std::min(gaussian_mult,
                     get_curvature_multiplicator(
                         F, working_edge.get_point_A(), e_size, avg_e_size));
      }
      if (!F.get_gradient_at_point(working_edge.get_point_B()).is_zero()) {
        gaussian_mult =
            std::min(gaussian_mult,
                     get_curvature_multiplicator(
                         F, working_edge.get_point_B(), e_size, avg_e_size));
      }
      if (!F.get_gradient_at_point(midpoint).is_zero()) {
        gaussian_mult =
            std::min(gaussian_mult, get_curvature_multiplicator(
                                        F, midpoint, e_size, avg_e_size));
      }
      assertm(gaussian_mult != 10000, "3 singular points!");
      height = std::max(std::min(gaussian_mult, 1.3 * avg_e_size),
                        0.7 * avg_e_size) *
               sqrt(numeric(3.0)) / 2.0;
      //+avg_e_size * sqrt(numeric(3)) / 2 * 0.5;
    } else {
      std::cout << "no value opt proj center" << endl;
      height = edge_height;
    }

    const Vector direction = height * dir;
    assertm(direction * normal_at_edge_midpoint < kEps * height,
            "Wrong direction!");
    assertm(direction * Vector(working_edge.get_point_A(),
                               working_edge.get_point_B()) <
                kEps * height,
            "Wrong direction!");
    Point P(center, direction);
    assertm(Vector(center, P).get_length() - height < kEps,
            "Wrong point to project!");

    assertm(!F.get_gradient_at_point(P).is_zero(), "Zero gradient!");
    Vector normal = F.get_gradient_at_point(P).unit();
    std::optional<Point> projected = project(P, normal, F, e_size);
    return projected;
  }
  // end section -- change the incident face plane to tangent plane

  height = basic_height;
  Vector direction = height * find_direction(working_edge, incident_face);
  assertm(direction * incident_face.get_normal() < kEps, "Wrong direction!");
  assertm(direction * Vector(working_edge.get_point_A(),
                             working_edge.get_point_B()) <
              kEps,
          "Wrong direction!");

  Point P(center, direction);
  assertm(Vector(center, P).get_length() - height < kEps,
          "Wrong point to project!");

  assertm(!F.get_gradient_at_point(P).is_zero(), "Zero gradient!");
  Vector normal = F.get_gradient_at_point(P).unit();
  std::optional<Point> projected = project(P, normal, F, e_size);
  return projected;
}
// returns true if overlap triangle is added to mesh else returns false
bool BasicAlgorithm::fix_overlap(const HalfEdge &working_edge,
                                 const MeshPoint &overlap_point,
                                 const bool Delaunay) {
  // if Delaunay constraint is satisfied add the triangle to
  // triangulation and end
  if (!check_conditions(working_edge, overlap_point.get_point(), Delaunay,
                        overlap_point.get_index())) {
    return false;
  }
  /*assertm(!my_mesh.has_incoming_prev(working_edge, overlap_point) &&
              !my_mesh.has_outgoing_next(working_edge, overlap_point),
          "Creating overlap from non overlap point!");*/
  create_triangle(working_edge, overlap_point.get_point(), false,
                  overlap_point.get_index());
  // cout << "overlap" << endl;
  return true;
}

bool BasicAlgorithm::fix_close_points(const HalfEdge &working_edge,
                                      const Point &projected,
                                      const bool Delaunay) {
  /*
std::cout << "in close points" << endl;
auto surrounding_points = find_close_points(projected, working_edge);
std::cout << my_mesh.is_in_mesh(Edge(projected, working_edge.get_point_B()))
<< " "
<< my_mesh.is_in_mesh(Edge(working_edge.get_point_A(), projected))
<< endl;
std::cout << surrounding_points.has_value() << endl;
*/
  if (auto surrounding_points = find_close_points(projected, working_edge);
      surrounding_points.has_value()) {
    // points closer to projected point than 0.4*e_size sorted from closest
    vector<MeshPoint> close_points = surrounding_points.value();
    for (MeshPoint close_point : close_points) {
      // todo rewrite to reflect delaunay state
      /*
      std::cout << "distance betweend close points: "
                << Vector(close_point.get_point(), projected).get_length()
                << endl;
                */
      if (!my_mesh.is_boundary_point(close_point)) continue;
      if (check_conditions(working_edge, close_point.get_point(), Delaunay,
                           close_point.get_index())) {
        // fill it it's basic triangle
        MeshPointIndex opposite_point =
            my_mesh.get_halfedge(working_edge.get_next()).get_B();
        if (my_mesh.has_incoming_prev(working_edge, close_point.get_point()) &&
            my_mesh.has_outgoing_next(working_edge, close_point.get_point()) &&
            opposite_point != close_point.get_index()) {
          if (basic_triangle(working_edge, close_point)) {
            return true;
          }
        } else if (fix_overlap(working_edge, close_point, true)) {
          return true;
        }
      }
    }
  }
  return false;
}

bool BasicAlgorithm::fix_proj(const HalfEdge &working_edge,
                              const Point &projected, const bool Delaunay) {
  // TODO(kuska) cut triangle in half
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
  // std::cout << working_edge.get_index() << endl;
  /*
  Edge BC(working_edge.get_point_B(), projected);
  Edge CA(projected, working_edge.get_point_A());
  Edge BA(working_edge.get_point_B(), working_edge.get_point_A());
  Edge AC(working_edge.get_point_A(), projected);
  Edge CB(projected, working_edge.get_point_B());
  */
  /*
  if (my_mesh.is_in_mesh(BC))
    std::cout << "BC" << endl;
  else if (my_mesh.is_in_mesh(CA))
    std::cout << "CA" << endl;
  else if (my_mesh.is_in_mesh(BA))
    std::cout << "BA" << endl;
  else if (my_mesh.is_in_mesh(AC))
    std::cout << "AC" << endl;
  else if (my_mesh.is_in_mesh(CB))
    std::cout << "CB" << endl;
  else
    cout << "OK" << endl;
    */
  if (!check_conditions(working_edge, projected, Delaunay)) {
    return false;
  }
  assertm(!my_mesh.has_incoming_prev(working_edge, projected) &&
              !my_mesh.has_outgoing_next(working_edge, projected),
          "Creating new triangle with prev or next!");
  create_triangle(working_edge, projected, true);
  return true;
}

// checks if neighbour triangles have less than pi/2 angle
bool BasicAlgorithm::_check_normals(const HalfEdge &working_edge,
                                    const Point &point) const {
  assertm(working_edge.get_point_A() != point &&
              working_edge.get_point_B() != point,
          "Invalid triangle!");
  const Triangle new_triangle(working_edge.get_point_B(),
                              working_edge.get_point_A(), point);
  auto [prev_face_index, next_face_index] =
      find_neighbor_faces(working_edge, point);
  if (prev_face_index != kInvalidFaceIndex) {
    const Face &prev_face = my_mesh.get_face(prev_face_index);
    if (!neighbor_triangles_normal_check(new_triangle,
                                         prev_face.get_triangle())) {
      // std::cout << "Prev triangles don't pass normal check!" << endl;
      return false;
    }
  }
  if (next_face_index != kInvalidFaceIndex) {
    const Face &next_face = my_mesh.get_face(next_face_index);
    if (!neighbor_triangles_normal_check(new_triangle,
                                         next_face.get_triangle())) {
      // std::cout << "Next triangles don't pass normal check!" << endl;
      return false;
    }
  }
  if (!neighbor_triangles_normal_check(
          my_mesh.get_face(working_edge.get_incident()).get_triangle(),
          new_triangle)) {
    // std::cout << "Neighbor triangles don't pass normal check!" << endl;
    return false;
  }
  return true;
}

bool BasicAlgorithm::_fix_prev_next(const HalfEdge &working_edge,
                                    const MeshPoint &point,
                                    const bool Delaunay) {
  // potentialy new triangle
  assertm(working_edge.get_point_A() != point.get_point() &&
              working_edge.get_point_B() != point.get_point(),
          "Invalid triangle!");
  Triangle new_triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                        point);

  if (!check_conditions(working_edge, point.get_point(), Delaunay,
                        point.get_index())) {
    return false;
  }
  // TODO(kuska) change the type to bool - new/not new
  create_triangle(working_edge, point.get_point(), false, point.get_index());
  return true;
}

bool BasicAlgorithm::fix_prev_next(const HalfEdge &working_edge,
                                   const bool Delaunay) {
  vector<MeshPoint> prev = my_mesh.get_prev(working_edge);
  vector<MeshPoint> next = my_mesh.get_next(working_edge);
  const MeshPoint other_point = my_mesh.get_meshpoint(
      my_mesh.get_halfedge(working_edge.get_next()).get_B());
  if (prev.size() == 1 && next.size() == 1 && prev[0] == next[0] &&
      basic_triangle(working_edge, prev[0]))
    return true;
  if (prev.size() == 1 &&  // prev[0] != other_point &&
      _fix_prev_next(working_edge, prev[0], Delaunay)) {
    // std::cout << "Fix prev!" << endl;
    return true;
  }
  if (next.size() == 1 &&  // next[0] != other_point &&
      _fix_prev_next(working_edge, next[0], Delaunay)) {
    // std::cout << "Fix next!" << endl;
    return true;
  }
  return false;
}

bool BasicAlgorithm::fix_breakers(const HalfEdge &working_edge,
                                  const Point &projected, const bool Delaunay) {
  assertm(working_edge.get_point_A() != projected &&
              working_edge.get_point_B() != projected,
          "Invalid triangle!");
  Triangle proj_T(working_edge.get_point_B(), working_edge.get_point_A(),
                  projected);

  // points that break Delaunay constraint
  vector<MeshPoint> breakers = my_mesh.get_breakers(proj_T);
  // std::cout << "breakers count: " << breakers.size() << endl;
  //  sort from closest to working_edge
  std::sort(breakers.begin(), breakers.end(),
            [this, &working_edge](const auto &i, const auto &j) {
              return my_mesh.line_point_dist(working_edge, i) <
                     my_mesh.line_point_dist(working_edge, j);
            });

  //  try create triangle with breakers
  for (MeshPoint point : breakers) {
    if (my_mesh.is_boundary_point(point) &&
        fix_overlap(working_edge, point, Delaunay)) {
      // cout << "here" << endl;
      return true;
    }
  }
  return false;
}

// one step of the algorithm
bool BasicAlgorithm::step(const HalfEdgeIndex &working_edge_index) {
  const HalfEdge working_edge = my_mesh.get_halfedge(working_edge_index);

  assertm(!is_active(working_edge_index),
          "Working edge found in active edges!");
  assertm(!is_checked(working_edge_index),
          "Working edge found in checked edges!");
  assertm(!is_bounding(working_edge_index),
          "Working edge found in bounding edges!");
  assertm(!is_border(working_edge_index),
          "Working edge found in border edges!");
  // std::cout << "Before basic_triangle" << endl;
  //  if there is a hole in triangle shape, fill it with triangle
  const MeshPointIndex basic_triangle_index = find_triangle_hole(working_edge);
  if (basic_triangle_index != -1 &&
      basic_triangle(working_edge,
                     my_mesh.get_meshpoint(basic_triangle_index))) {
    // cout << "fix basic triangle" << endl;
    return true;
  }
  // std::cout << "After basic_triangle" << endl;

  const std::optional<Point> projected = get_projected(working_edge);
  // unable to project
  if (!projected.has_value()) {
    push_edge_to_checked(working_edge_index);
    return false;
  }
  std::optional<Point> clipped = bounding_box.crop_to_box(
      working_edge.get_midpoint(), projected.value(), e_size, F);
  assertm(clipped.has_value(), "Unable to crop to box!");
  if (!bounding_box.is_inside(clipped.value()) &&
      !bounding_box.is_on(clipped.value())) {
    push_edge_to_checked(working_edge_index);
    return false;
  }
  if (fix_close_points(working_edge, clipped.value(), true)) {
    return true;
  }

  if (fix_breakers(working_edge, projected.value(), true)) {
    return true;
  }

  if (fix_proj(working_edge, clipped.value(), true)) {
    return true;
  }

  if (fix_prev_next(working_edge, true)) {
    return true;
  }

  assertm(!is_border(working_edge_index),
          "Working edge found in border edges!");
  push_edge_to_checked(working_edge_index);
  // std::cout << "Pushed edge to checked!" << endl;
  return false;
}

void BasicAlgorithm::starting() {
  int round = 0;
  int active_rounds = 0;
  while (!my_mesh.active_edges_empty()) {
    // my_mesh.obj_format(name);
    HalfEdgeIndex working_edge = kInvalidEdgeIndex;
    working_edge = round;  // my_mesh.get_active_edge();
    if (!my_mesh.is_active(working_edge)) {
      round++;
      continue;
    }
    assertm(my_mesh.is_in_mesh(my_mesh.get_halfedge(working_edge).get_edge()),
            "Working edge not in mesh edges!");
    assertm(my_mesh.is_boundary(working_edge), "Working edge not boundary!");
    assertm(my_mesh.get_halfedge(working_edge).is_active(),
            "Working edge not active!");

    assertm(!my_mesh.is_in_mesh(
                Edge(my_mesh.get_halfedge(working_edge).get_point_B(),
                     my_mesh.get_halfedge(working_edge).get_point_A())),
            "Opposite edge found in mesh esges!");
    my_mesh.remove_active_edge(working_edge);

    assertm(working_edge != kInvalidEdgeIndex, "Active edge is invalid!");
    //  step returns true if new triangle is created
    if (step(working_edge)) {
      assertm(!my_mesh.is_boundary(working_edge),
              "Working edge found in border!");

    } else {
      assertm(my_mesh.is_checked(working_edge), "Checked edge not checked");
    }

    // once in every 10 steps prints number of triangles and number of active
    // edges and updates output file
    if (active_rounds % 50 == 0) {
      my_mesh.cout_triangles_number();
      my_mesh.obj_format(name);

      cout << "Number of active edges: " << my_mesh.get_active_edges_size()
           << endl;
      cout << "Number of checked edges: " << my_mesh.get_checked_edges_size()
           << endl
           << endl;
    }
    round++;
    active_rounds++;
  }

  // output
  // my_mesh.obj_format(name);
  return;
}

void BasicAlgorithm::ending() {
  assertm(my_mesh.active_edges_empty(),
          "Called ending with non-empty active edges!");

  // else push checked edges to active, clear checked and call ending

  int round0 = 0;
  while (!my_mesh.checked_edges_empty() && round0 < 2) {
    round0++;
    my_mesh.move_active_edges_to_checked();

    int round = 0;
    // cout<<"Fixing holes!"<<endl;
    while (!my_mesh.active_edges_empty()) {
      round++;
      HalfEdgeIndex working_edge = kInvalidEdgeIndex;
      working_edge = my_mesh.get_active_edge();
      my_mesh.remove_active_edge(working_edge);

      fix_holes2(my_mesh.get_halfedge(working_edge));
      my_mesh.obj_format(name);
      if (round % 5 == 0) {
        my_mesh.cout_triangles_number();
        cout << "Number of active edges: " << my_mesh.get_active_edges_size()
             << endl
             << endl;
      }
    }
  }
  // my_mesh.obj_format(name);
  if (!my_mesh.checked_edges_empty()) {
    cout << "Some holes stayed!" << endl;
  }
  return;
}

// TODO(kuska) make this function void
void BasicAlgorithm::calculate() {
  starting();
  // my_mesh.obj_format(name + "_before_fix");
  std::cout << "Number of active after starting: "
            << my_mesh.get_active_edges_size() << endl;
  std::cout << "Number of checked after starting: "
            << my_mesh.get_checked_edges_size() << endl;
  ending();
  // my_mesh.obj_format(name + "_after_fix");
  std::cout << "Number of active after ending: "
            << my_mesh.get_active_edges_size() << endl;
  std::cout << "Number of checked after ending: "
            << my_mesh.get_checked_edges_size() << endl;
  // fix_corners();
  //   my_mesh.measure(bounding_edges, F, name, e_size);
  //   my_mesh.adaptive(0.005, F, e_size);
  //  my_mesh.mesh_format(name);
  cout << "Final triangles count: ";
  my_mesh.cout_triangles_number();
  my_mesh.obj_format(name);
  return;
}

// private member functions

// TODO(kuska) rewrite ending
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

// TODO(kuska) finish fix corners
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

// TODO(kuska) rewrite hole fixing algorithm

bool BasicAlgorithm::fix_holes2(const HalfEdge &working_edge) {
  // if there is a hole in triangle shape, fill it with triangle
  const MeshPointIndex basic_triangle_index = find_triangle_hole(working_edge);
  if (basic_triangle_index != -1 &&
      basic_triangle(working_edge,
                     my_mesh.get_meshpoint(basic_triangle_index))) {
    // cout << "fix basic triangle" << endl;
    return true;
  }

  //  find candidate point for working_edge
  std::optional<Point> projected = get_projected(working_edge);
  if (!projected.has_value()) {
    push_edge_to_checked(working_edge.get_index());
    return false;
  }
  std::optional<Point> clipped =
      bounding_box.crop_to_box(working_edge.get_midpoint(), projected.value(),
                               working_edge.get_length(), F);
  assertm(clipped.has_value(), "Unable to crop to box!");
  if (!bounding_box.is_inside(clipped.value()) &&
      !bounding_box.is_on(clipped.value())) {
    push_edge_to_checked(working_edge.get_index());
    return false;
  }
  /*
  std::optional<MeshPoint> closest_point =
      get_closest_point(working_edge);
  if (closest_point.has_value()) {
    create_triangle(working_edge, closest_point.value().get_point(), "prev",
                    closest_point.value().get_index());
    return true;
  }
  */
  if (fix_close_points(working_edge, clipped.value(), false)) {
    // cout << "fix close points" << endl;
    return true;
  }

  if (fix_breakers(working_edge, projected.value(), false)) {
    // cout << "fix breakers" << endl;
    return true;
  }
  /*x
  if (fix_proj(working_edge, projected, false)) {
    // cout << "fix proj" << endl;
    return true;
  }
*/
  if (fix_prev_next(working_edge, false)) {
    return true;
  }
  push_edge_to_checked(working_edge.get_index());
  return false;
}
