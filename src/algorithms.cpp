#include "algorithms.h"

#include <ginac/ginac.h>

#include <utility>

#include "assertm.h"
#include "function.h"
#include "point.h"

using namespace GiNaC;
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
    iter -= ex_to<numeric>(f.subs(my_x == iter).evalf() /
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
      0)
    projected = Bisect(my_x, f, new_point1, last_point1, 0);
  else if (f.subs(my_x == new_point2).evalf() *
               f.subs(my_x == last_point2).evalf() <=
           0)
    projected = Bisect(my_x, f, new_point2, last_point2, 0);
  else
    assertm(false, "Wrong call in Bisection function!");

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

  ex f = F.get_function().subs(
      lst{F.get_x() == param_x, F.get_y() == param_y, F.get_z() == param_z});
  ex df = f.diff(my_x, 1);

  std::optional<Point> projected = std::nullopt;

  numeric root = Newton_Raphson(my_x, f, df, starting_point);
  numeric projected_x = ex_to<numeric>(param_x.subs(my_x == root).evalf());
  numeric projected_y = ex_to<numeric>(param_y.subs(my_x == root).evalf());
  numeric projected_z = ex_to<numeric>(param_z.subs(my_x == root).evalf());

  projected = Point(projected_x, projected_y, projected_z);
  if (e_size.has_value() &&
      Vector(point_to_project, projected.value()).get_length() >
          4 * e_size.value()) {
    root = Bisection(my_x, f, starting_point, e_size.value());

    numeric projected_x = ex_to<numeric>(param_x.subs(my_x == root).evalf());
    numeric projected_y = ex_to<numeric>(param_y.subs(my_x == root).evalf());
    numeric projected_z = ex_to<numeric>(param_z.subs(my_x == root).evalf());
    assertm(F.substitute(lst{F.get_x() == projected_x, F.get_y() == projected_y,
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
/*vector<Edge> connect_edges(const vector<Edge> &v1, const vector<Edge> &v2) {
  vector<Edge> connected = v1;
  for (auto edge : v2) {
    connected.push_back(edge);
  }
  return connected;
}

// connects two vectors of points
vector<Point> connect_points(const vector<Point> &v1, const vector<Point> &v2) {
  vector<Point> connected = v1;
  for (auto point : v2) {
    connected.push_back(point);
  }
  return connected;
}
*/
// angle BAP in range (-Pi, Pi) with respect to neighbour triangle
numeric angle(const Edge &working_edge, const Point &P, const Triangle &N) {
  if (P == working_edge.A() || P == working_edge.B()) return false;
  Triangle T = Triangle(working_edge.A(), working_edge.B(), P);
  if (!T.is_triangle()) return false;
  assertm(N.is_triangle(), "Invalid neighbour triangle!");

  // third point in neighbour triangle

  std::optional<Point> NPoint = std::nullopt;

  if (N.A() != working_edge.A() && N.A() != working_edge.B())
    NPoint = N.A();

  else if (N.B() != working_edge.A() && N.B() != working_edge.B())
    NPoint = N.B();

  else if (N.C() != working_edge.A() && N.C() != working_edge.B())
    NPoint = N.C();

  assertm(NPoint.has_value(), "Not found neighbour point!");

  Vector AB = Vector(working_edge.A(), working_edge.B());
  Vector AP = Vector(working_edge.A(), P);
  Vector AN = Vector(working_edge.A(), NPoint.value());

  assertm(!AB.is_zero() && !AP.is_zero() && !AN.is_zero(), "Zero vectors!");

  numeric cos = (AB.unit() * AP.unit());
  assertm(abs(cos) <= numeric(1), "Wrong value of cosine!");

  numeric angle = acos(cos);

  assertm(angle >= 0, "Angle in the wrong range!");

  Vector ABcrossAN = (AB ^ AN);
  Vector ABcrossAP = (AB ^ AP);

  // if same_direction is > 0 the normals are pointing in aprox same direction
  // thus the points lie on the same side

  numeric same_direction = ABcrossAP * ABcrossAN;

  if (same_direction > 0) {
    angle = -angle;
  }
  assertm(abs(angle) <= ex_to<numeric>(Pi.evalf()),
          "Angle in the wrong range!");

  return angle;
}

// true if angle is between 0 and 3*pi/4 with respect to neighbour triangle
bool good_orientation(const Edge &working_edge, const Point &P,
                      const Triangle &N) {
  return angle(working_edge, P, N) > 0 &&
         angle(working_edge, P, N) < 3 * ex_to<numeric>(Pi.evalf()) / 4;
}

// https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d

// returns ditance between point and line segment given by working edge
numeric line_point_dist(const Edge &working_edge, const Point &P,
                        const Triangle &neighbour_triangle) {
  assertm(!Vector(working_edge.A(), working_edge.B()).is_zero(),
          "Edge is zero vector!");
  Vector AB_unit = Vector(working_edge.A(), working_edge.B()).unit();
  Vector AP = Vector(working_edge.A(), P);
  numeric t = AP * AB_unit;
  Point p = Point(working_edge.A(), t * AB_unit);

  numeric angle1 = angle(working_edge, P, neighbour_triangle);
  numeric angle2 =
      angle(Edge(working_edge.B(), working_edge.A()), P, neighbour_triangle);

  if (abs(angle1) <= ex_to<numeric>(Pi.evalf()) / 2 &&
      abs(angle2) <= ex_to<numeric>(Pi.evalf()) / 2) {
    return Vector(P, p).get_length();
  } else if (abs(angle1) > ex_to<numeric>(Pi.evalf()) / 2) {
    return Vector(working_edge.A(), P).get_length();
  } else {
    return Vector(working_edge.B(), P).get_length();
  }
  assertm(false, "Wrong line point distance!");
  return 1000;
}

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(const Edge &e, const Triangle &T) {
  assertm(T.is_triangle(), "Getting normal of non-valid triangle!");
  Vector normal = T.get_normal();
  Vector edge_vector(e.A(), e.B());
  assertm(!(normal ^ edge_vector).is_zero(), "Zero vector!");
  Vector direction = (normal ^ edge_vector).unit();
  numeric min_side_length = std::min(
      T.AB().get_length(), std::min(T.BC().get_length(), T.CA().get_length()));
  numeric min_median_length =
      std::min(Vector(T.A(), T.BC().get_midpoint()).get_length(),
               std::min(Vector(T.B(), T.CA().get_midpoint()).get_length(),
                        Vector(T.C(), T.AB().get_midpoint()).get_length()));
  numeric delta = min_median_length / 20;

  bool repeat = true;
  int round = 0;
  while (repeat && round < 10) {
    Point P1 = Point(e.get_midpoint(), delta * direction);
    delta = delta / 2;
    round++;
    if (T.is_in_triangle(P1)) {
      if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction))) {
        direction = numeric(-1) * direction;
        repeat = false;
      } else
        assertm(false, "Both points in triangle!");
    } else if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction))) {
      repeat = true;
    } else {
      repeat = false;
    }
  }
  assertm(!repeat, "No Points in triangle!");

  return direction;
}
