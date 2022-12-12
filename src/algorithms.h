#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <optional>
#include <vector>

#include "assertm.h"
#include "edge.h"
#include "function.h"
#include "mesh.h"
#include "point.h"
#include "triangle.h"

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

// returns projected point in the direction of normal
Point project(const Point &point_to_project, const Vector &normal,
              const Function &F, const std::optional<numeric> e_size);

// connects two vectors of edges
vector<HalfEdgeIndex> connect_edges(const vector<HalfEdgeIndex> &v1,
                                    const vector<HalfEdgeIndex> &v2);

// connects two vectors of points
// vector<Point> connect_points(const vector<Point> &v1, const vector<Point>
// &v2);

// angle BAP in range (-Pi, Pi) with respect to neighbour triangle
std::optional<numeric> angle(const Mesh &mesh, const HalfEdge &working_edge,
                             const Point &P, const Face &incident_face,
                             const bool clockwise);

// true if angle is between 0 and 3*pi/4 with respect to neighbour triangle
bool good_orientation(const Mesh &mesh, const HalfEdge &working_edge,
                      const Point &P, const Face &incident_face);

// https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
// returns ditance between point and line given by working edge
numeric line_point_dist(const Mesh &mesh, const HalfEdge &working_edge,
                        const Point &P, const Face &incident_face);

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(const HalfEdge &working_edge, const Face &F);

#endif
