#include "algorithms.h"

#include <ginac/ginac.h>

#include <utility>

using GiNaC::ex;
using GiNaC::numeric;
using GiNaC::realsymbol;
using std::pair;

// N-R method for root finding, not necessarily the closest root
numeric Newton_Raphson(const realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point) {
  numeric iter = starting_point;
  int iterations = 0;
  while (abs(f.subs(my_x == iter).evalf()) > kEps && iterations < 10) {
    iterations++;
    assertm(abs(df.subs(my_x == iter).evalf()) > kEps * kEps,
            "Division by 0 in N-R method!");
    iter -= GiNaC::ex_to<numeric>(f.subs(my_x == iter).evalf() /
                                  df.subs(my_x == iter).evalf());
  }
  return iter;
}

// bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter) {
  if (f.subs(my_x == point1).evalf() < kEps * kEps) return point1;
  if (f.subs(my_x == point2).evalf() < kEps * kEps) return point2;
  assertm(iter < 1000, "Too much iterations in bisection method!");
  assertm(f.subs(my_x == point1).evalf() * f.subs(my_x == point2).evalf() <= 0,
          "Wrong call for Bisect function!");
  numeric new_point = (point1 + point2) / 2;
  std::optional<numeric> projected = std::nullopt;
  if (f.subs(my_x == point1).evalf() * f.subs(my_x == new_point).evalf() <= 0) {
    projected = Bisect(my_x, f, point1, new_point, ++iter);
  } else if (f.subs(my_x == point2).evalf() *
                 f.subs(my_x == new_point).evalf() <=
             0) {
    projected = Bisect(my_x, f, new_point, point2, ++iter);
  } else {
    assertm(false, "Wrong value in Bisect!");
  }
  assertm(projected.has_value(), "Projected point withou value!");
  return projected.value();
}

// Bisection is called when N-R proejcts to distant point
// finds 2 points on opposite sides of surface and returns result of bisection
// on these two points
std::optional<numeric> Bisection(const realsymbol my_x, const ex &f,
                                 numeric starting_point, numeric e_size) {
  numeric dx = e_size / 10;
  numeric new_point1 = starting_point, last_point1 = starting_point;
  numeric new_point2 = starting_point, last_point2 = starting_point;
  int iterations = 0;
  while ((f.subs(my_x == new_point1).evalf() *
                  f.subs(my_x == last_point1).evalf() >
              0 &&
          f.subs(my_x == new_point2).evalf() *
                  f.subs(my_x == last_point2).evalf() >
              0) &&
         iterations < 1000) {
    iterations++;
    last_point1 = new_point1;
    new_point1 += dx;
    last_point2 = new_point2;
    new_point2 -= dx;
  }
  // error while bisection
  if (iterations >= 1000) return std::nullopt;

  std::optional<numeric> projected = std::nullopt;

  if (f.subs(my_x == new_point1).evalf() *
          f.subs(my_x == last_point1).evalf() <=
      0) {
    projected = Bisect(my_x, f, new_point1, last_point1, 0);
  } else if (f.subs(my_x == new_point2).evalf() *
                 f.subs(my_x == last_point2).evalf() <=
             0) {
    projected = Bisect(my_x, f, new_point2, last_point2, 0);
  } else {
    assertm(false, "Wrong call in Bisection function!");
  }

  assertm(projected.has_value(), "Projected point without value!");
  return projected.value();
}

// rotates point about the axis given by the point rot_center and the vector u 
// by an angle
Point rotate(const Point &rot_center, const Point &to_rotate, const Vector &u,
             const numeric angle) {
  // rotating matrix around ROT_direction
  // cos(pi/3)+ux^2*(1-cos(pi/3))          ux*uy*(1-cos(pi/3)-uz*sin(pi/3)
  // ux*uz(1-cos(pi/3))+uy*sin(pi/3) ux*uy*(1-cos(pi/3)+uz*sin(pi/3)
  // cos(pi/3)+uy^2*(1-cos(pi/3))         uy*uz(1-cos(pi/3))-ux*sin(pi/3)
  // ux*uz(1-cos(pi/3))-uy*sin(pi/3)       uy*uz(1-cos(pi/3))+ux*sin(pi/3)
  // cos(pi/3)+uz^2*(1-cos(pi/3))

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
  const Vector start_vector(rot_center, to_rotate);
  const Vector end_vector(start_vector * rot_vector_x,
                          start_vector * rot_vector_y,
                          start_vector * rot_vector_z);
  return Point(rot_center, end_vector);
}

// performs circular bisection starting at start_point in the range of end_angle
// the rotation is about the vector u
Point circular_bisection(const Function &F, const Point &singular_point,
                         const Point &start_point, numeric end_angle,
                         const Vector &u /*rotation direction*/,
                         numeric e_size) {
  const Point end_point = rotate(singular_point, start_point, u, end_angle);
  assertm(
      GiNaC::ex_to<numeric>(F.get_function()
                                .subs(GiNaC::lst{F.get_x() == start_point.x(),
                                                 F.get_y() == start_point.y(),
                                                 F.get_z() == start_point.z()})
                                .evalf()) *
              GiNaC::ex_to<numeric>(
                  F.get_function()
                      .subs(GiNaC::lst{F.get_x() == end_point.x(),
                                       F.get_y() == end_point.y(),
                                       F.get_z() == end_point.z()})
                      .evalf()) <=
          0,
      "Wrong call for circular bisect function!");

  std::pair<Point, Point> bisection_param = {start_point, end_point};
  numeric my_angle = end_angle;
  for (int i = 0; i < 100; i++) {
    bisection_param =
        _circular_bisect(F, bisection_param.first, bisection_param.second,
                         singular_point, u, my_angle);
    if (bisection_param.first == bisection_param.second) {
      return bisection_param.first;
    }
    my_angle = my_angle / 2;
  }
  assertm(false, "Too much iterations in bisection method!");
  return start_point;
}

// called from circilar_bisection, returns pair of points for the next iteration
std::pair<Point, Point> _circular_bisect(const Function &F,
                                         const Point &start_point,
                                         const Point &end_point,
                                         const Point &singular_point,
                                         const Vector &u, const numeric angle) {
  const numeric subs_start =
      GiNaC::ex_to<numeric>(F.get_function()
                                .subs(GiNaC::lst{F.get_x() == start_point.x(),
                                                 F.get_y() == start_point.y(),
                                                 F.get_z() == start_point.z()})
                                .evalf());
  const numeric subs_end =
      GiNaC::ex_to<numeric>(F.get_function()
                                .subs(GiNaC::lst{F.get_x() == end_point.x(),
                                                 F.get_y() == end_point.y(),
                                                 F.get_z() == end_point.z()})
                                .evalf());
  if (abs(subs_start) < kEps * kEps) {
    return {start_point, start_point};
  }
  if (abs(subs_end) < kEps * kEps) {
    return {end_point, end_point};
  }
  const Point &mid_point = rotate(singular_point, start_point, u, angle / 2);

  const numeric subs_mid =
      GiNaC::ex_to<numeric>(F.get_function()
                                .subs(GiNaC::lst{F.get_x() == mid_point.x(),
                                                 F.get_y() == mid_point.y(),
                                                 F.get_z() == mid_point.z()})
                                .evalf());
  if (abs(subs_mid) < kEps * kEps) {
    return {mid_point, mid_point};
  }
  if (subs_start * subs_mid < 0) {
    return {start_point, mid_point};
  }
  if (subs_mid * subs_end < 0) {
    return {mid_point, end_point};
  }
  assertm(false, "Wrong place!");
  return {start_point, end_point};
}

// returns projected point in the direction of normal
std::optional<Point> project(const Point &point_to_project,
                             const Vector &normal, const Function &F,
                             const numeric e_size /* = -1 */) {
  realsymbol my_x("my_x");

  numeric starting_point;

  Point P = point_to_project;
  Vector n = normal;

  ex param_x, param_y, param_z;

  // parametric equations of line given by P and n
  // expressing parameter and substituing to other equations

  if (abs(n.x()) > kEps) {
    starting_point = P.x();
    param_x = my_x;
    param_y = P.y() - n.y() * P.x() / n.x() + (n.y() / n.x()) * my_x;
    param_z = P.z() - n.z() * P.x() / n.x() + (n.z() / n.x()) * my_x;
  } else if (abs(n.y()) > kEps) {
    starting_point = P.y();
    param_x = P.x();
    param_y = my_x;
    param_z = P.z() - n.z() * P.y() / n.y() + (n.z() / n.y()) * my_x;
  } else if (abs(n.z()) > kEps) {
    starting_point = P.z();
    param_x = P.x();
    param_y = P.y();
    param_z = my_x;
  } else {
    assertm(false, "Normal is a zero vector!");
  }

  // after this in param_x, param_y and param_z are equations of only one
  // variable my_x substitute to F to get function for Newton-Raphson method

  ex f = F.get_function().subs(GiNaC::lst{
      F.get_x() == param_x, F.get_y() == param_y, F.get_z() == param_z});
  ex df = f.diff(my_x, 1);
  std::optional<Point> projected = std::nullopt;

  std::optional<numeric> root = Newton_Raphson(my_x, f, df, starting_point);
  assertm(root.has_value(), "Newton Rhapson not working!");
  numeric projected_x =
      GiNaC::ex_to<numeric>(param_x.subs(my_x == root.value()).evalf());
  numeric projected_y =
      GiNaC::ex_to<numeric>(param_y.subs(my_x == root.value()).evalf());
  numeric projected_z =
      GiNaC::ex_to<numeric>(param_z.subs(my_x == root.value()).evalf());

  projected = Point(projected_x, projected_y, projected_z);
  if (e_size != -1 &&
      Vector(point_to_project, projected.value()).get_length() > 4 * e_size) {
    root = Bisection(my_x, f, starting_point, e_size);
    if (!root.has_value()) return std::nullopt;

    numeric projected_x =
        GiNaC::ex_to<numeric>(param_x.subs(my_x == root.value()).evalf());
    numeric projected_y =
        GiNaC::ex_to<numeric>(param_y.subs(my_x == root.value()).evalf());
    numeric projected_z =
        GiNaC::ex_to<numeric>(param_z.subs(my_x == root.value()).evalf());
    assertm(F.substitute(GiNaC::lst{F.get_x() == projected_x,
                                    F.get_y() == projected_y,
                                    F.get_z() == projected_z}) < 10e-6,
            "Wrong Bisection calculation!");
    projected = Point(projected_x, projected_y, projected_z);
  }
  assertm(projected.has_value(), "Not found projected point!");
  return projected.value();
}

// finds the direction perpendicular to the working_edge in the plane of the face F 
// pointing outside of a face F
Vector find_direction_plane(const HalfEdge &working_edge, const Vector &normal,
                            const Face &F /*incident_face*/) {
  Vector edge_vector(working_edge.get_point_A(), working_edge.get_point_B());
  assertm(!(normal ^ edge_vector).is_zero(), "Zero vector!");
  Vector direction = (normal ^ edge_vector).unit();

  const Vector to_gravity_center(working_edge.get_midpoint(),
                                 F.get_gravity_center());
  if (direction * to_gravity_center > 0) direction = -1 * direction;

  return direction;
}

// returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(const HalfEdge &working_edge, const Face &F) {
  assertm(F.is_triangle(), "Getting normal of non-valid triangle!");
  Vector normal = F.get_normal();
  return find_direction_plane(working_edge, normal, F);
}

// curvedness of a surface F=0 in its point 
numeric get_curvedness(const Function &F, const Point &point) {
  vector<ex> gradient = F.get_gradient();
  Vector gradient_eval = F.get_gradient_at_point(point);
  assertm(gradient_eval != Vector(0, 0, 0), "zero gradient!");

  ex _Fxx = diff(gradient[0], F.get_x());
  ex _Fxy = diff(gradient[0], F.get_y());
  ex _Fxz = diff(gradient[0], F.get_z());
  ex _Fyx = diff(gradient[1], F.get_x());
  ex _Fyy = diff(gradient[1], F.get_y());
  ex _Fyz = diff(gradient[1], F.get_z());
  ex _Fzx = diff(gradient[2], F.get_x());
  ex _Fzy = diff(gradient[2], F.get_y());
  ex _Fzz = diff(gradient[2], F.get_z());

  // assertm(_Fxy == _Fyx && _Fyz == _Fzy && _Fxz == _Fzx, "Not smooth!");

  numeric Fxx = GiNaC::ex_to<numeric>(
      _Fxx.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());
  numeric Fxy = GiNaC::ex_to<numeric>(
      _Fxy.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());
  numeric Fxz = GiNaC::ex_to<numeric>(
      _Fxz.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());

  numeric Fyx = GiNaC::ex_to<numeric>(
      _Fyx.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());
  numeric Fyy = GiNaC::ex_to<numeric>(
      _Fyy.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());
  numeric Fyz = GiNaC::ex_to<numeric>(
      _Fyz.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());

  numeric Fzx = GiNaC::ex_to<numeric>(
      _Fzx.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());
  numeric Fzy = GiNaC::ex_to<numeric>(
      _Fzy.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());
  numeric Fzz = GiNaC::ex_to<numeric>(
      _Fzz.subs(GiNaC::lst{F.get_x() == point.x(), F.get_y() == point.y(),
                           F.get_z() == point.z()})
          .evalf());

  // hessian matrix
  Vector Gx(Fxx, Fxy, Fxz);
  Vector Gy(Fyx, Fyy, Fyz);
  Vector Gz(Fzx, Fzy, Fzz);

  // cofactor matrix
  Vector Hx(Fyy * Fzz - Fyz * Fzy, Fyz * Fzx - Fyx * Fzz,
            Fyx * Fzy - Fyy * Fzx);
  Vector Hy(Fxz * Fzy - Fxy * Fzz, Fxx * Fzz - Fxz * Fzx,
            Fxy * Fzx - Fxx * Fzy);
  Vector Hz(Fxy * Fyz - Fxz * Fyy, Fyx * Fxz - Fxx * Fzy,
            Fxx * Fyy - Fxy * Fyx);

  numeric gaussian_curvature =
      (gradient_eval *
       (Vector(gradient_eval * Hx, gradient_eval * Hy, gradient_eval * Hz))) /
      (gradient_eval.get_length_squared() * gradient_eval.get_length_squared());

  numeric mean_curvature =
      ((gradient_eval *
        (Vector(gradient_eval * Gx, gradient_eval * Gy, gradient_eval * Gz))) -
       gradient_eval.get_length_squared() * (Fxx + Fyy + Fzz)) /
      (2 * gradient_eval.get_length_squared() * gradient_eval.get_length());

  numeric sq;
  if (abs(mean_curvature * mean_curvature - gaussian_curvature) < kEps) {
    sq = 0;
  } else if (mean_curvature * mean_curvature - gaussian_curvature < -kEps) {
    assertm(false, "Incorrect curvatures!");
  } else {
    sq = sqrt(mean_curvature * mean_curvature - gaussian_curvature);
  }

  numeric k1, k2;
  if (abs(mean_curvature) < kEps && abs(gaussian_curvature) < kEps) {
    k2 = 0;
    k1 = 0;
  } else if (abs(gaussian_curvature) < kEps) {
    k2 = 0;
    k1 = mean_curvature * 2;
  } else {
    k1 = -mean_curvature + sq;
    k2 = -mean_curvature - sq;
  }
  numeric curvedness = sqrt((k1 * k1 + k2 * k2) / 2);
  return curvedness;
}

// returns the adaptive height
numeric get_curvature_multiplicator(
    const Function &F, const Point &point, const numeric &e_size,
    const numeric &working_edge_length) {
  // e_size is edge size for unit sphere
  // curvedness is the inverse of the radius of the sphere with given curvedness
  const numeric k =
      1.0 /
      e_size;  // how many times does the edge fit into radius of unit sphere
  numeric min_size = 0.7 * working_edge_length;
  numeric max_size = 1.3 * working_edge_length;

  const numeric curvedness = get_curvedness(F, point);
  if (curvedness < kEps) return max_size;
  numeric radius = 1.0 / curvedness;
  if (min_size > radius / k)
    return min_size;
  else if (max_size < radius / k)
    return max_size;
  return radius / k;
}

// finds first edge from seed point
Edge get_seed_edge(Point seed_point, const Function &F, numeric edge_size,
                   const BoundingBox &bounding_box) {
  numeric edge_length =
      get_curvature_multiplicator(F, seed_point, edge_size, edge_size);
  edge_length = edge_size; // comment this line for adaptive triangulation
  Vector edge_size_tangent =
      std::max(std::min(edge_length, edge_size * 1.3), edge_size * 0.7) *
      sqrt(3) / 2 * (F.get_tangent_at_point(seed_point).unit());

  Point point_to_project(seed_point, edge_size_tangent);

  // direction of projection
  Vector direction = F.get_gradient_at_point(point_to_project).unit();

  std::optional<Point> projected_point =
      project(point_to_project, direction, F, edge_length);
  assertm(projected_point.has_value(), "No value!");
  projected_point = bounding_box.crop_to_box(
      seed_point, projected_point.value(), edge_length, F);

  assertm(seed_point != projected_point.value(), "Error in get_seed_edge");

  return Edge(seed_point, projected_point.value());
}

// finds third point in first triangle from seed edge
Point get_seed_triangle(const Edge &e, numeric edge_size, const Function &F,
                        const BoundingBox &bounding_box) {
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
  numeric basic_height = sqrt(numeric(3)) / 2 * e.get_length();

  const Point mult_point = Point(basic_height * center_tangent.unit().x(),
                                 basic_height * center_tangent.unit().y(),
                                 basic_height * center_tangent.unit().z());
  numeric edge_length_1 =
      get_curvature_multiplicator(F, e.A(), edge_size, e.get_length());
  numeric edge_length_2 =
      get_curvature_multiplicator(F, e.B(), edge_size, e.get_length());
  numeric edge_length_3 = get_curvature_multiplicator(
      F, mult_point, edge_size, e.get_length());
  numeric edge_length =
      std::min(std::min(edge_length_1, edge_length_2), edge_length_3);
  edge_length = edge_size; // comment this line for adaptive triangulation
  numeric height =
      sqrt(numeric(3)) / 2 *
      std::max(std::min(edge_length, edge_size * 1.3), edge_size * 0.7);

  Point point_to_project(center, height * center_tangent.unit());

  Vector normal = F.get_gradient_at_point(point_to_project).unit();

  std::optional<Point> projected =
      project(point_to_project, normal, F, edge_length);
  assertm(projected.has_value(), "No value!");
  projected = bounding_box.crop_to_box(e.get_midpoint(), projected.value(),
                                       edge_length, F);
  return projected.value();
}

// returns first triangle
Triangle find_seed_triangle(const Function &F, Point seed, numeric e_size,
                            const BoundingBox &bounding_box) {
  Vector normal = F.get_gradient_at_point(seed).unit();
  // project point on surface just to be sure it is lying on the surface with
  // enough precision
  std::optional<Point> proj_seed = project(seed, normal, F, e_size);
  assertm(proj_seed.has_value(), "No value!");
  seed = proj_seed.value();
  assertm(bounding_box.is_inside(seed), "Seed point outside of bounding box!");
  // gets seed edge
  Edge seed_edge = get_seed_edge(seed, F, e_size, bounding_box);

  // gets third point in seed triangle
  Point Q = get_seed_triangle(seed_edge, e_size, F, bounding_box);

  // return seed triangle
  return Triangle(seed_edge.A(), seed_edge.B(), Q);
}