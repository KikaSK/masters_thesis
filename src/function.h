#pragma once

#include <ginac/ginac.h>

#include <iostream>
#include <vector>

#include "assertm.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"

using GiNaC::ex;
using GiNaC::numeric;
using GiNaC::realsymbol;
using std::endl;
using std::vector;

class Function {
 private:
  const realsymbol _x;
  const realsymbol _y;
  const realsymbol _z;
  const ex _F;
  const vector<ex> _dF;

 public:
  Function(const realsymbol &x, const realsymbol &y, const realsymbol &z, ex F,
           vector<ex> dF);
  Function() = default;

  ex get_function() const;
  vector<ex> get_gradient() const;
  realsymbol get_x() const;
  realsymbol get_y() const;
  realsymbol get_z() const;

  ex grad_x() const;
  ex grad_y() const;
  ex grad_z() const;

  Vector get_gradient_at_point(const Point &P) const;
  Vector get_tangent_at_point(const Point &P) const;
  numeric eval_at_point(const Point &P) const;
  bool is_inside(const Point &P) const;
  bool is_on(const Point &P) const;
  bool is_outside(const Point &P) const;
  Vector outside_normal(const Triangle &T, const numeric e_size) const;
  numeric substitute(GiNaC::ex il) const;
};
