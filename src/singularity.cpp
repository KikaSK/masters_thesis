#include "singularity.h"

Singularity::Singularity(const Point location, const vector<Vector>& directions,
                         const SingularityType type, const int n)
    : _location(location), _directions(directions), _type(type), _n(n){};

numeric Singularity::x() const { return _location.x(); }
numeric Singularity::y() const { return _location.y(); }
numeric Singularity::z() const { return _location.z(); }
Point Singularity::location() const { return _location; }
vector<Vector> Singularity::directions() const { return _directions; }
Vector Singularity::get_direction(int i) const {
  assertm(i < _directions.size(), "Invalid index!");
  return _directions[i];
}
int Singularity::get_directions_count() const { return _directions.size(); }
SingularityType Singularity::type() const { return _type; }
int Singularity::n() const { return _n; }
void Singularity::set_location(Point location) { _location = location; }
void Singularity::set_directions(vector<Vector> directions) {
  _directions = directions;
}
void Singularity::add_direction(Vector direction) {
  _directions.push_back(direction);
}
void Singularity::set_type(SingularityType type) { _type = type; }
void Singularity::set_n(int n) { _n = n; }
