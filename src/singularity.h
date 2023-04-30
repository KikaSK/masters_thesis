#pragma once

#include <ginac/ginac.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "assertm.h"
#include "constants.h"
#include "vector.h"

using GiNaC::numeric;
using std::endl;
using std::vector;

enum class SingularityType {
  Amm = 0,  // A--
  Apm = 1,  // A+-
  Dmm = 2,  // D--
  Dpm = 3,  // D+-
  Epp = 4,  // E++
  Epm = 5,  // E+-
};

class Singularity {
 private:
  Point _location;
  vector<Vector> _directions;
  SingularityType _type;
  int _n;

  // TODO Ajko rule of 3

 public:
  Singularity(const Point location, const vector<Vector>& directions,
              const SingularityType type, const int n);
  Singularity() = default;

  // getters
  numeric x() const;
  numeric y() const;
  numeric z() const;
  Point location() const;
  vector<Vector> directions() const;
  Vector get_direction(int i) const;
  int get_directions_count() const;
  SingularityType type() const;
  int n() const;

  // setters
  void set_location(Point location);
  void set_directions(vector<Vector> directions);
  void add_direction(Vector direction);
  void set_type(SingularityType type);
  void set_n(int n);
};
