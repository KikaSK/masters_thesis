#pragma once

#include <optional>
#include <vector>

#include "assertm.h"
#include "bounding_box.h"
#include "edge.h"
#include "face.h"
#include "function.h"
#include "half_edge.h"
#include "point.h"
#include "triangle.h"

using GiNaC::ex;
using GiNaC::numeric;
using GiNaC::realsymbol;
using std::vector;

class BoundingBox;

// N-R method for root finding, not necessarily the closest root
numeric Newton_Raphson(const realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point);

// bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter);

// Bisection is called when N-R proejcts to distant point
// finds 2 points on opposite sides of surface and returns result of bisection
// on these two points
std::optional<numeric> Bisection(const realsymbol my_x, const ex &f,
                                 numeric starting_point, numeric e_size);

// performs circular bisection starting at start_point in the range of end_angle
// the rotation is about the vector u
Point circular_bisection(const Function &F, const Point &singular_point,
                         const Point &start_point, numeric end_angle,
                         const Vector &u /*rotation direction*/,
                         numeric e_size);

// called from circilar_bisection, returns pair of points for the next iteration
std::pair<Point, Point> _circular_bisect(const Function &F,
                                         const Point &start_point,
                                         const Point &end_point,
                                         const Point &singular_point,
                                         const Vector &u, const numeric angle);

// rotates point about the axis given by the point rot_center and the vector u 
// by an angle
Point rotate(const Point &rot_center, const Point &to_rotate, const Vector &u,
             const numeric angle);

// returns projected point in the direction of normal
std::optional<Point> project(const Point &point_to_project,
                             const Vector &normal, const Function &F,
                             const numeric e_size = -1);

// finds the direction perpendicular to the working_edge in the plane of the face F 
// pointing outside of a face F
Vector find_direction_plane(const HalfEdge &working_edge, const Vector &normal,
                            const Face &F);

// returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(const HalfEdge &working_edge, const Face &F);

// returns the adaptive height
numeric get_curvature_multiplicator(
    const Function &F, const Point &point, const numeric &e_size,
    const numeric &working_edge_length);

// returns first triangle
Triangle find_seed_triangle(const Function &F, Point seed, numeric e_size,
                            const BoundingBox &bounding_box);
