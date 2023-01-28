#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <optional>
#include <vector>

#include "assertm.h"
#include "edge.h"
#include "function.h"
//#include "mesh.h"
#include "face.h"
#include "half_edge.h"
#include "point.h"
#include "triangle.h"

using GiNaC::ex;
using GiNaC::numeric;
using GiNaC::realsymbol;
using std::vector;

// N-R method for root finding, not necessarily the closest root
numeric Newton_Raphson(const realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point);

// bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter);

// Bisection is called when N-R proejcts to distant point
// finds 2 points on opposite sides of surface and returns result of bisection
// on these two points
numeric Bisection(const realsymbol my_x, const ex &f, numeric starting_point,
                  numeric e_size);

Point circular_bisection(const Function &F, const Point &singular_point,
                         const Point &start_point, numeric end_angle,
                         const Vector &u /*rotation direction*/,
                         numeric e_size);

std::pair<Point, Point> _circular_bisect(const Function &F,
                                         const Point &start_point,
                                         const Point &end_point,
                                         const Point &singular_point,
                                         const Vector &u, const numeric angle);

Point rotate(const Point &rot_center, const Point &to_rotate, const Vector &u,
             const numeric angle);

// returns projected point in the direction of normal
Point project(const Point &point_to_project, const Vector &normal,
              const Function &F, const std::optional<numeric> e_size);

// connects two vectors of edges
vector<HalfEdgeIndex> connect_edges(const vector<HalfEdgeIndex> &v1,
                                    const vector<HalfEdgeIndex> &v2);

// connects two vectors of points
// vector<Point> connect_points(const vector<Point> &v1, const vector<Point>
// &v2);

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(const HalfEdge &working_edge, const Face &F);

numeric get_gaussian_curvature_multiplicator(const Function &F,
                                             const Point &point);

#endif
