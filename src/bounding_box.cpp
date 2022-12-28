#include "bounding_box.h"

BoundingBox::BoundingBox(numeric _min_x, numeric _max_x, numeric _min_y,
                         numeric _max_y, numeric _min_z, numeric _max_z)
    : _min_x(_min_x),
      _max_x(_max_x),
      _min_y(_min_y),
      _max_y(_max_y),
      _min_z(_min_z),
      _max_z(_max_z){};

numeric BoundingBox::min_x() const { return _min_x; }
numeric BoundingBox::max_x() const { return _max_x; }
numeric BoundingBox::min_y() const { return _min_y; }
numeric BoundingBox::max_y() const { return _max_y; }
numeric BoundingBox::min_z() const { return _min_z; }
numeric BoundingBox::max_z() const { return _max_z; }

bool BoundingBox::in_interval_x(const numeric x) const {
  return (x >= _min_x && x <= _max_x);
}
bool BoundingBox::in_interval_y(const numeric y) const {
  return (y >= _min_y && y <= _max_y);
}
bool BoundingBox::in_interval_z(const numeric z) const {
  return (z >= _min_z && z <= _max_z);
}

bool BoundingBox::is_inside(const Point P) const {
  return in_interval_x(P.x()) && in_interval_y(P.y()) && in_interval_z(P.z());
}

bool BoundingBox::is_on(const Point P) const {
  bool x_wall =
      ((abs(P.x() - _min_x) < 10e-10 || abs(P.x() - _max_x) < 10e-10) &&
       in_interval_y(P.y()) && in_interval_z(P.z()));
  bool y_wall =
      ((abs(P.y() - _min_y) < 10e-10 || abs(P.y() - _max_y) < 10e-10) &&
       in_interval_x(P.x()) && in_interval_z(P.z()));
  bool z_wall =
      ((abs(P.z() - _min_z) < 10e-10 || abs(P.z() - _max_z) < 10e-10) &&
       in_interval_y(P.y()) && in_interval_x(P.x()));
  return (x_wall || y_wall || z_wall);
}

bool BoundingBox::new_bounding_edge(const Edge &e) const {
  return is_on(e.A()) && is_on(e.B());
}

int BoundingBox::faces(const Point &P) const {
  int output = 0;
  if (abs(P.x() - _min_x) < 10e-5) output |= 1;
  if (abs(P.x() - _max_x) < 10e-5) output |= 1 << 1;
  if (abs(P.y() - _min_y) < 10e-5) output |= 1 << 2;
  if (abs(P.y() - _max_y) < 10e-5) output |= 1 << 3;
  if (abs(P.z() - _min_z) < 10e-5) output |= 1 << 4;
  if (abs(P.z() - _max_z) < 10e-5) output |= 1 << 5;
  return output;
}

std::optional<Point> BoundingBox::project_on_min_x(const Point &midpoint,
                                                   const Point &P) const {
  Vector v = Vector(P, midpoint);
  // Vector v = Vector(1, 0, 0);

  numeric dist = Vector(midpoint, P).get_length();

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_x - _min_x;
  vector<ex> d_input = {1, 0, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.x() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) ||
        Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_max_x(const Point &midpoint,
                                                   const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  // Vector v = Vector(1, 0, 0);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_x - _max_x;
  vector<ex> d_input = {1, 0, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.x() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) ||
        Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_min_y(const Point &midpoint,
                                                   const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  // Vector v = (0, 1, 0);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_y - _min_y;
  vector<ex> d_input = {0, 1, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.y() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) ||
        Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_max_y(const Point &midpoint,
                                                   const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  // Vector v = (0, 1, 0);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_y - _max_y;
  vector<ex> d_input = {0, 1, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.y() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) ||
        Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_min_z(const Point &midpoint,
                                                   const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  // Vector v = Vector(0, 0, 1);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_z - _min_z;
  vector<ex> d_input = {0, 0, 1};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.z() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) ||
        Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_max_z(const Point &midpoint,
                                                   const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  // Vector v = Vector(0, 0, 1);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_z - _max_z;
  vector<ex> d_input = {0, 0, 1};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.z() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) ||
        Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}

std::optional<Point> BoundingBox::project_on_box(const Point &midpoint,
                                                 const Point &P) const {
  // projecting point on all of the walls
  vector<std::optional<Point>> points = {std::nullopt, std::nullopt,
                                         std::nullopt, std::nullopt,
                                         std::nullopt, std::nullopt};
  points[0] = project_on_min_x(midpoint, P);
  points[1] = project_on_max_x(midpoint, P);
  points[2] = project_on_min_y(midpoint, P);
  points[3] = project_on_max_y(midpoint, P);
  points[4] = project_on_min_z(midpoint, P);
  points[5] = project_on_max_z(midpoint, P);

  std::optional<numeric> min_dist = std::nullopt;
  std::optional<int> min_dist_index = std::nullopt;
  for (int i = 0; i < 6; ++i) {
    if (points[i].has_value()) {
      numeric dist = Vector(P, points[i].value()).get_length();
      if (!min_dist.has_value()) {
        min_dist = dist;
        min_dist_index = i;
      } else if (dist < min_dist.value()) {
        min_dist = dist;
        min_dist_index = i;
      }
    }
  }
  // assertm(min_dist_index.has_value(), "Index without value!");
  if (!min_dist_index.has_value()) return std::nullopt;
  assertm(points[min_dist_index.value()].has_value(),
          "Projected point without value!");
  return points[min_dist_index.value()].value();
}

std::optional<Point> BoundingBox::crop_to_box(const Point &midpoint,
                                              const Point &P,
                                              const numeric &e_size,
                                              const Function &F) const {
  numeric precision = e_size / 4;
  std::optional<Point> projected = P;

  if (!is_inside(P)) {
    projected = project_on_box(midpoint, P);
  } else if (abs(P.x() - _min_x) < precision &&
             project_on_min_x(midpoint, P).has_value()) {
    projected = project_on_min_x(midpoint, P).value();
  } else if (abs(P.x() - _max_x) < precision &&
             project_on_max_x(midpoint, P).has_value()) {
    projected = project_on_max_x(midpoint, P).value();
  } else if (abs(P.y() - _min_y) < precision &&
             project_on_min_y(midpoint, P).has_value()) {
    projected = project_on_min_y(midpoint, P).value();
  } else if (abs(P.y() - _max_y) < precision &&
             project_on_max_y(midpoint, P).has_value()) {
    projected = project_on_max_y(midpoint, P).value();
  } else if (abs(P.z() - _min_z) < precision &&
             project_on_min_z(midpoint, P).has_value()) {
    projected = project_on_min_z(midpoint, P).value();
  } else if (abs(P.z() - _max_z) < precision &&
             project_on_max_z(midpoint, P).has_value()) {
    projected = project_on_max_z(midpoint, P).value();
  } else {
    projected = P;
  }
  if (!projected.has_value()) return std::nullopt;
  int faces_index = faces(projected.value());
  Vector v(1, 1, 1);
  Point line_p(0, 0, 0);
  bool close_wall = false;

  // lies on min_x wall
  if (faces_index & 1) {
    v = v - Vector(1, 0, 0);
    line_p = Point(_min_x, line_p.y(), line_p.z());
    if (abs(projected.value().y() - _min_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _min_y, line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().y() - _max_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _max_y, line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().z() - _min_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _min_z);
      close_wall = true;
    }
    if (abs(projected.value().z() - _max_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _max_z);
      close_wall = true;
    }
  }
  // lies on max_x wall
  if (faces_index & 1 << 1) {
    v = v - Vector(1, 0, 0);
    line_p = Point(_max_x, line_p.y(), line_p.z());
    if (abs(projected.value().y() - _min_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _min_y, line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().y() - _max_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _max_y, line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().z() - _min_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _min_z);
      close_wall = true;
    }
    if (abs(projected.value().z() - _max_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _max_z);
      close_wall = true;
    }
  }
  // lies on min_y wall
  if (faces_index & 1 << 2) {
    v = v - Vector(0, 1, 0);
    line_p = Point(line_p.x(), _min_y, line_p.z());
    if (abs(projected.value().x() - _min_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_min_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().x() - _max_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_max_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().z() - _min_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _min_z);
      close_wall = true;
    }
    if (abs(projected.value().z() - _max_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _max_z);
      close_wall = true;
    }
  }
  // lies on max_y wall
  if (faces_index & 1 << 3) {
    v = v - Vector(0, 1, 0);
    line_p = Point(line_p.x(), _max_y, line_p.z());
    if (abs(projected.value().x() - _min_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_min_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().x() - _max_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_max_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().z() - _min_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _min_z);
      close_wall = true;
    }
    if (abs(projected.value().z() - _max_z) < precision) {
      v = v - Vector(0, 0, 1);
      line_p = Point(line_p.x(), line_p.y(), _max_z);
      close_wall = true;
    }
  }
  // lies on min_z wall
  if (faces_index & 1 << 4) {
    v = v - Vector(0, 0, 1);
    line_p = Point(line_p.x(), line_p.y(), _min_z);
    if (abs(projected.value().x() - _min_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_min_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().x() - _max_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_max_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().y() - _min_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _min_y, line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().y() - _max_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _max_y, line_p.z());
      close_wall = true;
    }
  }
  // lies on max_z wall
  if (faces_index & 1 << 5) {
    v = v - Vector(0, 0, 1);
    line_p = Point(line_p.x(), line_p.y(), _max_z);
    if (abs(projected.value().x() - _min_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_min_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().x() - _max_x) < precision) {
      v = v - Vector(1, 0, 0);
      line_p = Point(_max_x, line_p.y(), line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().y() - _min_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _min_y, line_p.z());
      close_wall = true;
    }
    if (abs(projected.value().y() - _max_y) < precision) {
      v = v - Vector(0, 1, 0);
      line_p = Point(line_p.x(), _max_y, line_p.z());
      close_wall = true;
    }
  }
  if (close_wall) {
    if (line_p.x() == 0)
      line_p = Point(projected.value().x(), line_p.y(), line_p.z());
    if (line_p.y() == 0)
      line_p = Point(line_p.x(), projected.value().y(), line_p.z());
    if (line_p.z() == 0)
      line_p = Point(line_p.x(), line_p.y(), projected.value().z());

    std::optional<Point> clipped_point = std::nullopt;

    if (v.is_zero()) {
      clipped_point = projected.value();
    } else {
      clipped_point = project(line_p, v, F, e_size);
    }
    assertm(clipped_point.has_value(), "Point without value!");
    if (Vector(clipped_point.value(), projected.value()).get_length() <
        e_size) {
      return clipped_point.value();
    }
  }
  return projected.value();
}
