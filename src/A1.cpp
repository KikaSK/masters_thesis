#pragma once

#include "basic_algorithm.h"

/*
void BasicAlgorithm::triangulate_A1_starter(const Function &F, const numeric
&e_size, const Point &singular, const Vector &direction) { Vector unit_direction
= direction.unit(); Vector plane_vector =
unit_direction.get_any_perpendicular(); Point point_to_project = Point(singular,
(sqrt(15)/4)*plane_vector + (e_size/4)*unit_direction);
    // direction of projection
    Vector direction = F.get_gradient_at_point(point_to_project).unit();
    Point projected_point = project(point_to_project, direction, F, {e_size});

    // checks if the point is in the correct halfspace
    Vector projected_vector = Vector(singular, projected_point);
    assertm(direction*projected_vector > 0, "New point in opposite halfspace!");

    const Edge new_edge_1(singular, projected_point);

    Point center = new_edge_1.get_midpoint();

    // normals at endpoints of the edge
    Vector center_normal = F.get_gradient_at_point(center).unit();

    const Vector &edge_vector = projected_vector;
    Vector center_tangent = center_normal ^ edge_vector;

    assertm(abs(center_normal * center_tangent) < 1e-6, "Not perpendicular!");

    // height of equilateral triangle with side edge_size
    numeric height = sqrt(numeric(3)) / 2 * e_size;
    Point point_to_project_2(center, height * center_tangent.unit());
    Vector normal = F.get_gradient_at_point(point_to_project_2).unit();
    Point projected = project(point_to_project_2, normal, F, {e_size});
}*/