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
  numeric precision = 1e-6;

  numeric iter = starting_point;
  int iterations = 0;
  while (abs(f.subs(my_x == iter).evalf()) > precision && iterations < 1'000) {
    iterations++;
    assertm(abs(df.subs(my_x == iter).evalf()) > 10e-6,
            "Division by 0 in N-R method!");
    iter -= GiNaC::ex_to<numeric>(f.subs(my_x == iter).evalf() /
                                  df.subs(my_x == iter).evalf());
  }
  return iter;
}

// bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter) {
  if (f.subs(my_x == point1).evalf() < 10e-6) return point1;
  if (f.subs(my_x == point2).evalf() < 10e-6) return point2;
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
numeric Bisection(const realsymbol my_x, const ex &f, numeric starting_point,
                  numeric e_size) {
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

  assertm(iterations < 1000,
          "Bisection method failed, try smaller triangle edge size");

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

// returns projected point in the direction of normal
Point project(const Point &point_to_project, const Vector &normal,
              const Function &F, const std::optional<numeric> e_size) {
  realsymbol my_x("my_x");

  numeric starting_point;

  Point P = point_to_project;
  Vector n = normal;

  ex param_x, param_y, param_z;

  // parametric equations of line given by P and n
  // expressing parameter and substituing to other equations

  if (n.x() != 0) {
    starting_point = P.x();
    param_x = my_x;
    param_y = P.y() - n.y() * P.x() / n.x() + (n.y() / n.x()) * my_x;
    param_z = P.z() - n.z() * P.x() / n.x() + (n.z() / n.x()) * my_x;
  } else if (n.y() != 0) {
    starting_point = P.y();
    param_x = P.x();
    param_y = my_x;
    param_z = P.z() - n.z() * P.y() / n.y() + (n.z() / n.y()) * my_x;
  } else if (n.z() != 0) {
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

  numeric root = Newton_Raphson(my_x, f, df, starting_point);
  numeric projected_x =
      GiNaC::ex_to<numeric>(param_x.subs(my_x == root).evalf());
  numeric projected_y =
      GiNaC::ex_to<numeric>(param_y.subs(my_x == root).evalf());
  numeric projected_z =
      GiNaC::ex_to<numeric>(param_z.subs(my_x == root).evalf());

  projected = Point(projected_x, projected_y, projected_z);
  if (e_size.has_value() &&
      Vector(point_to_project, projected.value()).get_length() >
          4 * e_size.value()) {
    root = Bisection(my_x, f, starting_point, e_size.value());

    numeric projected_x =
        GiNaC::ex_to<numeric>(param_x.subs(my_x == root).evalf());
    numeric projected_y =
        GiNaC::ex_to<numeric>(param_y.subs(my_x == root).evalf());
    numeric projected_z =
        GiNaC::ex_to<numeric>(param_z.subs(my_x == root).evalf());
    assertm(F.substitute(GiNaC::lst{F.get_x() == projected_x,
                                    F.get_y() == projected_y,
                                    F.get_z() == projected_z}) < 10e-6,
            "Wrong Bisection calculation!");
  }
  assertm(projected.has_value(), "Not found projected point!");
  /*
  assertm(!e_size.has_value() ||
              Vector(point_to_project, projected.value()).get_length() <
                  4 * e_size.value(),
          "Wrong calculation in project function!");
  */
  return projected.value();
}

// connects two vectors of edges
vector<HalfEdgeIndex> connect_edges(const vector<HalfEdgeIndex> &v1,
                                    const vector<HalfEdgeIndex> &v2) {
  vector<HalfEdgeIndex> connected = v1;
  for (auto edge : v2) {
    connected.push_back(edge);
  }
  return connected;
}
/*
// connects two vectors of points
vector<Point> connect_points(const vector<Point> &v1, const vector<Point> &v2) {
  vector<Point> connected = v1;
  for (auto point : v2) {
    connected.push_back(point);
  }
  return connected;
}
*/

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(const HalfEdge &working_edge, const Face &F) {
  assertm(F.is_triangle(), "Getting normal of non-valid triangle!");
  Vector normal = F.get_normal();
  Vector edge_vector(working_edge.get_point_A(), working_edge.get_point_B());
  assertm(!(normal ^ edge_vector).is_zero(), "Zero vector!");
  Vector direction = (normal ^ edge_vector).unit();
  numeric min_side_length = std::min(
      F.AB().get_length(), std::min(F.BC().get_length(), F.CA().get_length()));
  numeric min_median_length =
      std::min(Vector(F.A(), F.BC().get_midpoint()).get_length(),
               std::min(Vector(F.B(), F.CA().get_midpoint()).get_length(),
                        Vector(F.C(), F.AB().get_midpoint()).get_length()));
  numeric delta = min_median_length / 20;

  bool repeat = true;
  int round = 0;
  while (repeat && round < 10) {
    Point P1 = Point(working_edge.get_midpoint(), delta * direction);
    delta = delta / 2;
    round++;
    if (F.is_in_triangle(P1)) {
      if (!F.is_in_triangle(
              Point(working_edge.get_midpoint(), -delta * direction))) {
        direction = numeric(-1) * direction;
        repeat = false;
      } else {
        assertm(false, "Both points in triangle!");
      }
    } else if (!F.is_in_triangle(
                   Point(working_edge.get_midpoint(), -delta * direction))) {
      repeat = true;
    } else {
      repeat = false;
    }
  }
  assertm(!repeat, "No Points in triangle!");
  return direction;
}
