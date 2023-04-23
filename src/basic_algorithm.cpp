#include "basic_algorithm.h"

#include <iostream>

using std::cout;

#pragma region "Region singularities"

void BasicAlgorithm::get_local_mesh(Mesh *local_mesh, const Function &_F,
                                    const Function &_G,
                                    const Function &inter_FG,
                                    const vector<Point> &polyline,
                                    const BoundingBox &bounding_box,
                                    const numeric e_size) {
  assertm(polyline.size() > 1, "Wrong length of polyline!");

  const Function *F = &_F;
  const Function *G = &_G;

  int offset_index_points = 0;
  int offset_index_edges = 0;

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
    if (F->is_outside(tan_G_point)) tangent_G = -1 * tangent_G;
    if (G->is_outside(tan_F_point)) tangent_F = -1 * tangent_F;
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

    if (test_T.get_normal() * inter_FG.outside_normal(test_T, e_size) < 0) {
      std::swap(projected_F, projected_G);
      std::swap(F, G);
    }

    const Triangle T_F(A, B, projected_F);
    const Triangle T_G(B, A, projected_G);
    if (i == 1) {
      local_mesh->add_first_triangle(T_F, bounding_box);
    } else if (i == polyline.size() - 1 && polyline[0] == polyline[i]) {
      // special case - closed curve
      local_mesh->add_triangle_to_two_meshpoints(
          offset_index_points + 1 + (i - 2) * 3, 0, projected_F, bounding_box);
    } else {
      local_mesh->add_triangle_to_meshpoint(
          offset_index_points + 1 + (i - 2) * 3, B, projected_F, bounding_box);
    }

    local_mesh->add_triangle(offset_index_edges + (i - 1) * 6, projected_G,
                             true, bounding_box);
  }
  local_mesh->obj_format("non-isolated-test");
  return;
}

void BasicAlgorithm::triangulate_An_analytical(
    const Point &singular, const Vector &singular_direction, Mesh *mesh,
    const MeshPointIndex singular_index, const int n) {
  Vector unit_singular_direction = singular_direction.unit();
  Vector plane_vector = unit_singular_direction.get_any_perpendicular().unit();
  const Vector &u = unit_singular_direction;
  const numeric e = e_size;
  const numeric h = pow(e_size, (numeric)2 / (n + 1));
  int num_triangles = 6;
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

  vector<Point> points;
  Vector rot_vector = (plane_vector ^ unit_singular_direction).unit();
  Vector direction = plane_vector * e + u * h;
  for (int i = 0; i < num_triangles; ++i) {
    Point new_point =
        Point(singular.x() + direction.x(), singular.y() + direction.y(),
              singular.z() + direction.z());
    points.push_back(new_point);
    // checks if the point is in the correct halfspace
    assertm(unit_singular_direction * direction >= 0,
            "New point in opposite halfspace!");
    direction = Vector(rot_vector_x * direction, rot_vector_y * direction,
                       rot_vector_z * direction);
  }
  const Triangle first_triangle(singular, points[0], points[1]);

  // make sure normals point outside
  if (first_triangle.get_normal() * F.outside_normal(first_triangle, e_size) <
      0)
    reverse(points.begin(), points.end());

  // Mesh local_mesh;
  for (int i = 0; i < num_triangles; ++i) {
    int j = (i + 1) % num_triangles;
    Triangle triangle = Triangle(singular, points[i], points[j]);
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (mesh->get_faces_count() == 0 ||
          singular_index == kInvalidPointIndex) {
        mesh->add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        mesh->add_triangle_to_meshpoint(singular_index, points[i], points[j],
                                        bounding_box);
      }
    } else if (i < num_triangles - 1) {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], true,
                         bounding_box);
      mesh->edges_check("in A1 round: " + std::to_string(i));
      //,bounding_box);
    } else {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], false,
                         bounding_box,
                         mesh->get_points_count() - num_triangles);
      mesh->edges_check("in A1 round: " + std::to_string(i));
    }
  }

  mesh->obj_format(name + "_local");
  return;
}

// triangulates singularities lying in a single half-space
void BasicAlgorithm::triangulate_singularity_circular(
    const Point &singular, const Vector &singular_direction, Mesh *mesh,
    const MeshPointIndex singular_index) {
  Vector unit_singular_direction = singular_direction.unit();
  Vector plane_vector = unit_singular_direction.get_any_perpendicular().unit();
  const Vector &u = unit_singular_direction;

  // rotating matrix around singular direction
  // cos(pi/3)+ux^2*(1-cos(pi/3))          ux*uy*(1-cos(pi/3)-uz*sin(pi/3)
  // ux*uz(1-cos(pi/3))+uy*sin(pi/3) ux*uy*(1-cos(pi/3)+uz*sin(pi/3)
  // cos(pi/3)+uy^2*(1-cos(pi/3))         uy*uz(1-cos(pi/3))-ux*sin(pi/3)
  // ux*uz(1-cos(pi/3))-uy*sin(pi/3)       uy*uz(1-cos(pi/3))+ux*sin(pi/3)
  // cos(pi/3)+uz^2*(1-cos(pi/3))

  const int num_triangles = 6;

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

  const numeric mult = 0.8;
  vector<Point> points;
  for (int i = 0; i < num_triangles; ++i) {
    Point start_point_to_bisect(singular, e_size * mult * plane_vector);
    Point end_point_to_bisect(singular,
                              e_size * mult * unit_singular_direction);
    Vector rot_vector = (plane_vector ^ unit_singular_direction).unit();
    Point projected_point =
        circular_bisection(F, singular, start_point_to_bisect,
                           GiNaC::ex_to<numeric>((GiNaC::Pi / 2).evalf()),
                           rot_vector, e_size * mult);
    points.push_back(projected_point);
    // checks if the point is in the correct halfspace
    Vector projected_vector = Vector(singular, projected_point);
    assertm(unit_singular_direction * projected_vector >= 0,
            "New point in opposite halfspace!");

    plane_vector =
        Vector(rot_vector_x * plane_vector, rot_vector_y * plane_vector,
               rot_vector_z * plane_vector);
  }
  assertm(
      singular != points[0] && singular != points[1] && points[0] != points[1],
      "Wrong projection!");
  const Triangle first_triangle(singular, points[0], points[1]);
  assertm(first_triangle.is_triangle(), "Not a triangle!");

  // make sure normals point outside - outside normal problem in higher degrees
  if (first_triangle.get_normal() * F.outside_normal(first_triangle, e_size) <
      0)
    reverse(points.begin(), points.end());

  // Mesh local_mesh;
  for (int i = 0; i < num_triangles; ++i) {
    int j = (i + 1) % num_triangles;
    assertm(singular != points[i] && singular != points[j] &&
                points[i] != points[j],
            "Wrong projection!");
    Triangle triangle = Triangle(singular, points[i], points[j]);
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (mesh->get_faces_count() == 0 ||
          singular_index == kInvalidPointIndex) {
        mesh->add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        mesh->add_triangle_to_meshpoint(singular_index, points[i], points[j],
                                        bounding_box);
      }
    } else if (i < num_triangles - 1) {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], true,
                         bounding_box);
      mesh->edges_check("in A1 round: " + std::to_string(i));
      //,bounding_box);
    } else {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], false,
                         bounding_box,
                         mesh->get_points_count() - num_triangles);
      mesh->edges_check("in A1 round: " + std::to_string(i));
    }
  }

  mesh->obj_format(name + "_local");
  return;
}

// trinagulates singularities lying in both half-spaces having only single
// branch
void BasicAlgorithm::triangulate_singularity_case2(
    const Point &singular, const Vector &singular_direction, Mesh *mesh,
    const MeshPointIndex singular_index) {
  Vector unit_singular_direction = singular_direction.unit();
  Vector plane_vector = unit_singular_direction.get_any_perpendicular().unit();
  const Vector &u = unit_singular_direction;
  const int num_triangles = 4;
  std::cout << "ok1" << endl;
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
  std::cout << "ok2" << endl;
  const numeric mult = 0.7;
  vector<Point> points;
  int extra_triangles = 0;
  for (int i = 0; i < num_triangles + extra_triangles; ++i) {
    Point start_point_to_bisect(singular,
                                e_size * mult * unit_singular_direction);
    Point end_point_to_bisect(singular,
                              -e_size * mult * unit_singular_direction);
    Vector rot_vector = (plane_vector).unit();
    Point projected_point = circular_bisection(
        F, singular, start_point_to_bisect,
        GiNaC::ex_to<numeric>((GiNaC::Pi).evalf()), rot_vector, e_size * mult);
    points.push_back(projected_point);
    std::cout << projected_point << endl;

    // checks if the point is in the correct halfspace
    // Vector projected_vector = Vector(singular, projected_point);
    // assertm(unit_singular_direction * projected_vector > 0,
    //        "New point in opposite halfspace!");

    plane_vector =
        Vector(rot_vector_x * plane_vector, rot_vector_y * plane_vector,
               rot_vector_z * plane_vector);
  }

  std::cout << "ok3" << endl;
  assertm(
      singular != points[0] && singular != points[1] && points[0] != points[1],
      "Wrong projection!");
  const Triangle first_triangle(singular, points[0], points[1]);
  // make sure normals point outside
  if (first_triangle.get_normal() * F.outside_normal(first_triangle, e_size) <
      0)
    reverse(points.begin(), points.end());

  std::cout << "ok4" << endl;
  // Mesh local_mesh;
  for (int i = 0; i < num_triangles; ++i) {
    int j = (i + 1) % num_triangles;
    Triangle triangle = Triangle(singular, points[i], points[j]);
    assertm(singular != points[i] && singular != points[j] &&
                points[i] != points[j],
            "Wrong projection!");
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (mesh->get_faces_count() == 0) {
        mesh->add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        mesh->add_triangle_to_meshpoint(singular_index, points[i], points[j],
                                        bounding_box);
      }
    } else if (i < num_triangles - 1) {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], true,
                         bounding_box);
      mesh->edges_check("in A1 round: " + std::to_string(i));
      //,bounding_box);
    } else {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], false,
                         bounding_box,
                         mesh->get_points_count() - num_triangles);
      mesh->edges_check("in A1 round: " + std::to_string(i));
    }
  }
  std::cout << "ok4" << endl;
  mesh->obj_format(name + "_local");
  return;
}

void BasicAlgorithm::triangulate_cone_iterative(
    const Point &singular, const Vector &singular_direction, Mesh *mesh,
    const MeshPointIndex singular_index) {
  Vector unit_singular_direction = singular_direction.unit();
  Vector plane_vector = unit_singular_direction.get_any_perpendicular().unit();
  const Vector &u = unit_singular_direction;
  const int num_triangles = 4;
  const numeric mult = 1;
  vector<Point> points;
  const numeric delta_angle = GiNaC::ex_to<numeric>(GiNaC::Pi.evalf() / 36);
  for (int i = 0; i < num_triangles; ++i) {
    Point start_point_to_bisect(singular,
                                e_size * mult * unit_singular_direction);
    Point end_point_to_bisect = start_point_to_bisect;
    int iter = 0;
    while (F.eval_at_point(start_point_to_bisect) *
               F.eval_at_point(end_point_to_bisect) >
           0) {
      assertm(iter <= 71, "Too many iterations!");
      end_point_to_bisect =
          rotate(singular, end_point_to_bisect, plane_vector, delta_angle);
      iter++;
    }
    Point projected_point =
        circular_bisection(F, singular, start_point_to_bisect,
                           iter * delta_angle, plane_vector, e_size * mult);
    points.push_back(projected_point);
    plane_vector = Vector(
        singular,
        rotate(singular,
               Point(plane_vector.x() - singular.x(),
                     plane_vector.y() - singular.y(),
                     plane_vector.z() - singular.z()),
               u,
               GiNaC::ex_to<numeric>(2 * GiNaC::Pi.evalf() / num_triangles)));
  }

  const Triangle first_triangle(singular, points[0], points[1]);
  // make sure normals point outside
  if (first_triangle.get_normal() * F.outside_normal(first_triangle, e_size) <
      0)
    reverse(points.begin(), points.end());

  // Mesh local_mesh;
  for (int i = 0; i < num_triangles; ++i) {
    int j = (i + 1) % num_triangles;
    Triangle triangle = Triangle(singular, points[i], points[j]);
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (mesh->get_faces_count() == 0) {
        mesh->add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        mesh->add_triangle_to_meshpoint(singular_index, points[i], points[j],
                                        bounding_box);
      }
    } else if (i < num_triangles - 1) {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], true,
                         bounding_box);
      mesh->edges_check("in A1 round: " + std::to_string(i));
      //,bounding_box);
    } else {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], false,
                         bounding_box,
                         mesh->get_points_count() - num_triangles);
      mesh->edges_check("in A1 round: " + std::to_string(i));
    }
  }

  mesh->obj_format(name + "_local");
  return;
}

// TODO(kuska) project via spinning, not perpendicular
void BasicAlgorithm::triangulate_A1_starter(
    const Point &singular, const Vector &singular_direction, Mesh *mesh,
    const MeshPointIndex singular_index) {
  Vector unit_singular_direction = singular_direction.unit();
  Vector plane_vector = unit_singular_direction.get_any_perpendicular();
  const Vector &u = unit_singular_direction;

  // rotating matrix around singular direction
  // cos(pi/3)+ux^2*(1-cos(pi/3))          ux*uy*(1-cos(pi/3)-uz*sin(pi/3)
  // ux*uz(1-cos(pi/3))+uy*sin(pi/3) ux*uy*(1-cos(pi/3)+uz*sin(pi/3)
  // cos(pi/3)+uy^2*(1-cos(pi/3))         uy*uz(1-cos(pi/3))-ux*sin(pi/3)
  // ux*uz(1-cos(pi/3))-uy*sin(pi/3)       uy*uz(1-cos(pi/3))+ux*sin(pi/3)
  // cos(pi/3)+uz^2*(1-cos(pi/3))

  const int num_triangles = 5;

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

  vector<Point> points;
  for (int i = 0; i < num_triangles; ++i) {
    Point point_to_project =
        Point(singular, (sqrt(15) / 4) * plane_vector +
                            (e_size / 4) * unit_singular_direction);

    // direction of projection
    Vector direction = F.get_gradient_at_point(point_to_project).unit();
    std::optional<Point> projected_point =
        project(point_to_project, direction, F, e_size);
    assertm(projected_point.has_value(), "No value!");
    points.push_back(projected_point.value());

    // checks if the point is in the correct halfspace
    Vector projected_vector = Vector(singular, projected_point.value());
    assertm(direction * projected_vector > 0,
            "New point in opposite halfspace!");

    plane_vector =
        Vector(rot_vector_x * plane_vector, rot_vector_y * plane_vector,
               rot_vector_z * plane_vector);
  }

  const Triangle first_triangle(singular, points[0], points[1]);
  // make sure normals point outside
  if (first_triangle.get_normal() * F.outside_normal(first_triangle, e_size) <
      0)
    reverse(points.begin(), points.end());

  // Mesh local_mesh;
  for (int i = 0; i < num_triangles; ++i) {
    int j = (i + 1) % num_triangles;
    Triangle triangle = Triangle(singular, points[i], points[j]);
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (mesh->get_faces_count() == 0) {
        mesh->add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        mesh->add_triangle_to_meshpoint(singular_index, points[i], points[j],
                                        bounding_box);
      }
    } else if (i < num_triangles - 1) {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], true,
                         bounding_box);
      mesh->edges_check("in A1 round: " + std::to_string(i));
      //,bounding_box);
    } else {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], false,
                         bounding_box,
                         mesh->get_points_count() - num_triangles);
      mesh->edges_check("in A1 round: " + std::to_string(i));
    }
  }

  /*
  const numeric COS = numeric(1) / 2;
  const numeric SIN = numeric(sqrt(3) / 2);
  const Vector rot_vector_x = Vector(COS + u.x() * u.x() * (1 - COS),
                                     u.x() * u.y() * (1 - COS) - u.z() * SIN,
                                     u.x() * u.z() * (1 - COS) + u.y() * SIN);
  const Vector rot_vector_y = Vector(u.x() * u.y() * (1 - COS) + u.z() * SIN,
                                     COS + u.y() * u.y() * (1 - COS),
                                     u.y() * u.z() * (1 - COS) - u.x() * SIN);
  const Vector rot_vector_z = Vector(u.x() * u.z() * (1 - COS) - u.y() * SIN,
                                     u.y() * u.z() * (1 - COS) + u.x() * SIN,
                                     COS + u.z() * u.z() * (1 - COS));

  vector<Point> points;
  for (int i = 0; i < 6; ++i) {
    Point point_to_project =
        Point(singular, (sqrt(15) / 4) * plane_vector +
                            (e_size / 4) * unit_singular_direction);

    // direction of projection
    Vector direction = F.get_gradient_at_point(point_to_project).unit();
    Point projected_point = project(point_to_project, direction, F, {e_size});
    points.push_back(projected_point);

    // checks if the point is in the correct halfspace
    Vector projected_vector = Vector(singular, projected_point);
    assertm(direction * projected_vector > 0,
            "New point in opposite halfspace!");

    plane_vector =
        Vector(rot_vector_x * plane_vector, rot_vector_y * plane_vector,
               rot_vector_z * plane_vector);
  }

  // Mesh local_mesh;
  for (int i = 0; i < 6; ++i) {
    int j = (i + 1) % 6;
    Triangle triangle = Triangle(singular, points[i], points[j]);
    assertm(triangle.is_triangle(), "Invalid triangle!");
    if (i == 0) {
      if (mesh->get_faces_count() == 0) {
        mesh->add_first_triangle(triangle, bounding_box);
      } else {
        assertm(singular_index != kInvalidPointIndex,
                "invalid singular meshpoint index");
        mesh->add_triangle_to_meshpoint(singular_index, points[i], points[j],
                                        bounding_box);
      }
    } else if (i < 5) {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], true,
                         bounding_box);
      mesh->edges_check("in A1 round: " + std::to_string(i));
      //,bounding_box);
    } else {
      mesh->add_triangle(mesh->get_edges_count() - 1, points[j], false,
                         bounding_box, mesh->get_points_count() - 6);
      mesh->edges_check("in A1 round: " + std::to_string(i));
    }
  }
  */
  mesh->obj_format(name + "_local");
  return;
}

#pragma endregion "Region singularities"

#pragma region "Seed triangle"
// finds first edge from seed point
Edge BasicAlgorithm::get_seed_edge(Point seed_point) const {
  Vector edge_size_tangent =
      e_size * (F.get_tangent_at_point(seed_point).unit());

  Point point_to_project(seed_point, edge_size_tangent);

  // direction of projection
  Vector direction = F.get_gradient_at_point(point_to_project).unit();

  std::optional<Point> projected_point =
      project(point_to_project, direction, F, e_size);
  // projected_point = bounding_box.crop_to_box(seed_point, projected_point,
  // edge_size, F);
  assertm(projected_point.has_value(), "No value!");
  assertm(seed_point != projected_point.value(), "Error in get_seed_edge");

  return Edge(seed_point, projected_point.value());
}

// finds third point in first triangle from seed edge
Point BasicAlgorithm::get_seed_triangle(const Edge &e) const {
  Point center = e.get_midpoint();

  // normals at endpoints of the edge
  Vector n_A = F.get_gradient_at_point(e.A()).unit();
  Vector n_B = F.get_gradient_at_point(e.B()).unit();

  // average of endpoints normals
  Vector center_normal((n_A + n_B) / 2);

  Vector edge_vector(e.A(), e.B());
  Vector center_tangent = center_normal ^ edge_vector;

  assertm(abs(center_normal * center_tangent) < 1e-6, "Not perpendicular!");

  // height of equilateral triangle with side edge_size
  numeric height = sqrt(numeric(3)) / 2 * e_size;

  Point point_to_project(center, height * center_tangent.unit());

  Vector normal = F.get_gradient_at_point(point_to_project).unit();

  std::optional<Point> projected = project(point_to_project, normal, F, e_size);
  assertm(projected.has_value(), "No value!");
  // projected = bounding_box.crop_to_box(e.get_midpoint(), projected,
  // edge_size, F);
  return projected.value();
}

// returns first triangle
Triangle BasicAlgorithm::find_seed_triangle(Point seed) const {
  Vector normal = F.get_gradient_at_point(seed).unit();
  // project point on surface just to be sure it is lying on the surface with
  // enough precision
  std::optional<Point> opt_seed = project(seed, normal, F, e_size);
  assertm(opt_seed.has_value(), "No value!");
  seed = opt_seed.value();
  assertm(bounding_box.is_inside(seed), "Seed point outside of bounding box!");
  // gets seed edge
  Edge seed_edge = get_seed_edge(seed);

  // gets third point in seed triangle
  Point Q = get_seed_triangle(seed_edge);

  const Triangle new_triangle(seed_edge.A(), seed_edge.B(), Q);
  if (new_triangle.get_normal() * F.outside_normal(new_triangle, e_size) < 0)
    return Triangle(seed_edge.B(), seed_edge.A(), Q);

  // return seed triangle
  return new_triangle;
}
#pragma endregion "Seed triangle"

#pragma region "Find close points"

// TODO(kuska) rewrite to non-optional
std::optional<vector<MeshPoint>> BasicAlgorithm::_tree_close_points_finder(
    const Point &P, const HalfEdge &working_edge) const {
  numeric dist = 0.4 * e_size;
  auto close_points = my_mesh.get_meshpoints_in_interval(
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

// TODO(Kuska) add condition to check for the edge already existing
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
  // std::cout << "Check Delaunay conditions" << endl;
  if (working_edge.get_point_A() == P || working_edge.get_point_B() == P) {
    // std::cout << "Not triangle failed" << endl;
    return false;
  }
  Triangle new_triangle =
      Triangle(working_edge.get_point_B(), working_edge.get_point_A(), P);

  if (!new_triangle.is_triangle()) {
    // std::cout << "Is triangle failed" << endl;
    return false;
  }

  numeric min_edge_length = std::min(
      std::min(new_triangle.AB().get_length(), new_triangle.BC().get_length()),
      new_triangle.CA().get_length());

  if (my_mesh.get_faces_count() > 40 &&
      (new_triangle.AB().get_length() > 3 * min_edge_length ||
       new_triangle.CA().get_length() > 3 * min_edge_length ||
       new_triangle.BC().get_length() > 3 * min_edge_length)) {
    // std::cout << "Edge length failed" << endl;
    return false;
  }

  const Point opposite_point =
      my_mesh
          .get_meshpoint(my_mesh.get_halfedge(working_edge.get_next()).get_B())
          .get_point();

  // const Point midpoint = working_edge.get_midpoint();
  // const Vector edge_vector_new(midpoint, P);
  // const Vector edge_vector_neighbour(midpoint, opposite_point);

  if (P == opposite_point) return false;

  bool gc_distance = true;  // gc_distance_check(new_triangle);
  bool delaunay = my_mesh.check_Delaunay(working_edge, P, e_size);
  bool has_good_edges = good_edges(working_edge, P);
  bool neighbor_triangles_normal_check = _check_normals(working_edge, P);
  bool is_edge_in_mesh = my_mesh.is_in_mesh(new_triangle.AB()) ||
                         my_mesh.is_in_mesh(new_triangle.BC()) ||
                         my_mesh.is_in_mesh(new_triangle.CA());
  // TODO(kuska) Do not pass point to oreintability check!!!
  // bool orientability = orientability_check(working_edge, P);
  /*
  if (neighbor_triangles_normal_check) {
    cout << "Neighbor triangles normal check passed for triangles: " << endl;
    cout << new_triangle << endl;
    cout << "and" << endl;
    cout << my_mesh.get_face(working_edge.get_incident()).get_triangle()
         << endl;
  } else {
  */
  /*
    if (is_edge_in_mesh) cout << "NO PASS IS EDGE IN MESH!" << endl;
    if (!neighbor_triangles_normal_check) cout << "NO PASS NORMAL CHECK!" <<
    endl; if (!delaunay) cout << "NO PASS DELAUNAY!" << endl;
    // if (!orientability) cout << "NO PASS ORIENTABILITY!" << endl;
    if (!has_good_edges) cout << "NO PASS GOOD EDGES!" << endl;
  */
  return (delaunay && has_good_edges && neighbor_triangles_normal_check &&
          !is_edge_in_mesh && gc_distance
          // && orientability
  );
}

// checks if essential conditions are fulfilled
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
  //  if (!has_good_edges) cout << "Good Edges condition failed!" << endl;

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

  // std::cout << "in normal check:" << endl;

  Vector normal_1 = T1.get_normal();
  Vector normal_2 = T2.get_normal();

  // std::cout << "normal 1: " << normal_1 << endl;
  // std::cout << "normal 2: " << normal_2 << endl;
  // std::cout << "dot: " << normal_1 * normal_2 << endl;
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

bool BasicAlgorithm::overlap_normals_check(const MeshPoint &candidate,
                                           const HalfEdge &working_edge) const {
  if (candidate.get_point() == working_edge.get_point_A() ||
      candidate.get_point() == working_edge.get_point_B()) {
    return false;
  }

  Triangle my_triangle(working_edge.get_point_B(), working_edge.get_point_A(),
                       candidate.get_point());

  if (my_triangle.is_triangle() /*&&
      good_edges(working_edge, candidate.get_point())*/) {
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
  assertm(opposite_point != point, "opposite point");
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
  numeric basic_height = e_size * sqrt(numeric(3)) / 2;

  std::optional<Point> opt_proj_center =
      project(center, F.get_gradient_at_point(center).unit(), F, e_size);
  assertm(opt_proj_center.has_value(), "No value!");
  numeric gaussian_height =
      basic_height / get_curvature_multiplicator(F, opt_proj_center.value());

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
  numeric height = basic_height;
  //    gaussian_height * 0.5 + working_edge.get_length() * sqrt(numeric(3)) /
  //    4;
  //  height = gaussian_height;
  //   non adaptive height
  //   numeric height = gaussian_height;

  // section -- change the incident face plane to tangent plane

  const Face &incident_face = my_mesh.get_face(working_edge.get_incident());
  const Vector gradient_at_edge_midpoint = F.get_gradient_at_point(center);
  // Vector dir1(0, 0, 0);
  if (!gradient_at_edge_midpoint.is_zero()) {
    Vector normal_at_edge_midpoint = gradient_at_edge_midpoint.unit();
    if (normal_at_edge_midpoint * incident_face.get_normal() < 0)
      normal_at_edge_midpoint = -1 * normal_at_edge_midpoint;
    const Vector direction =
        height * find_direction_plane(working_edge, normal_at_edge_midpoint,
                                      incident_face);
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
    for (auto close_point : close_points) {
      // todo rewrite to reflect delaunay state
      /*
      std::cout << "distance betweend close points: "
                << Vector(close_point.get_point(), projected).get_length()
                << endl;
                */
      if (check_conditions(working_edge, close_point.get_point(), Delaunay,
                           close_point.get_index())) {
        // fill it it's basic triangle
        MeshPointIndex opposite_point =
            my_mesh.get_halfedge(working_edge.get_next()).get_B();
        if (my_mesh.has_incoming_prev(working_edge, close_point.get_point()) &&
            my_mesh.has_outgoing_next(working_edge, close_point.get_point()) &&
            opposite_point != close_point.get_index()) {
          if (basic_triangle(working_edge, close_point)) {
            // std::cout << "basic" << endl;
            return true;
          }
        } else if (fix_overlap(working_edge, close_point, true)) {
          // std::cout << "Does the point have good prev/next?" << endl;
          const MeshPoint &B = my_mesh.get_meshpoint(working_edge.get_B());
          const MeshPoint &A = my_mesh.get_meshpoint(working_edge.get_A());
          for (HalfEdgeIndex halfedge_index : B.get_outgoing()) {
            const HalfEdge &outgoing = my_mesh.get_halfedge(halfedge_index);
            if (outgoing.get_edge() ==
                Edge(B.get_point(), close_point.get_point())) {
              // cout << "outgoing next found!!!" << endl;
              break;
            }
          }
          for (HalfEdgeIndex halfedge_index : A.get_outgoing()) {
            const HalfEdgeIndex incoming_index =
                my_mesh.get_halfedge(halfedge_index).get_previous();
            const HalfEdge incoming = my_mesh.get_halfedge(incoming_index);
            if (incoming.get_edge() ==
                Edge(close_point.get_point(), A.get_point())) {
              // cout << "incoming prev found!!!" << endl;
              break;
            }
          }
          // std::cout << "overlap" << endl;
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
  for (auto point : breakers) {
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
    // cout << "working edge " << working_edge << endl;
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
    if (active_rounds % 30 == 0) {
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
  my_mesh.obj_format(name);

  /*
  // this means there are no holes so we finished
  if (checked_edges.empty()) {
    cout << "No holes!" << endl;
    return;
  }
  return;
  */
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
  if (!my_mesh.checked_edges_empty()) {
    cout << "Some holes stayed!" << endl;
  }
  return;
}

// TODO(kuska) make this function void
void BasicAlgorithm::calculate() {
  starting();
  my_mesh.obj_format(name + "_before_fix");
  std::cout << "Number of active after starting: "
            << my_mesh.get_active_edges_size() << endl;
  std::cout << "Number of checked after starting: "
            << my_mesh.get_checked_edges_size() << endl;
  ending();
  my_mesh.obj_format(name + "_after_fix");
  std::cout << "Number of active after ending: "
            << my_mesh.get_active_edges_size() << endl;
  std::cout << "Number of checked after ending: "
            << my_mesh.get_checked_edges_size() << endl;
  // fix_corners();
  //  my_mesh.measure(bounding_edges, F, name, e_size);
  //  my_mesh.adaptive(0.005, F, e_size);
  // my_mesh.mesh_format(name);
  cout << "Final triangles count: ";
  my_mesh.cout_triangles_number();
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
  std::optional<Point> clipped = bounding_box.crop_to_box(
      working_edge.get_midpoint(), projected.value(), e_size, F);
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
