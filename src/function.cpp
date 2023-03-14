#include "function.h"

Function::Function(const GiNaC::realsymbol &x, const GiNaC::realsymbol &y,
                   const GiNaC::realsymbol &z, ex F, vector<ex> dF)
    : _x(x), _y(y), _z(z), _F(F), _dF(dF){};

ex Function::get_function() const { return _F; }
vector<ex> Function::get_gradient() const { return _dF; }
realsymbol Function::get_x() const { return _x; }
realsymbol Function::get_y() const { return _y; }
realsymbol Function::get_z() const { return _z; }

ex Function::grad_x() const { return _dF[0]; }
ex Function::grad_y() const { return _dF[1]; }
ex Function::grad_z() const { return _dF[2]; }

Vector Function::get_gradient_at_point(const Point &P) const {
  // substituing point P to gradient
  numeric first = GiNaC::ex_to<numeric>(
      _dF[0].subs(GiNaC::lst{_x == P.x(), _y == P.y(), _z == P.z()}).evalf());
  numeric second = GiNaC::ex_to<numeric>(
      _dF[1].subs(GiNaC::lst{_x == P.x(), _y == P.y(), _z == P.z()}).evalf());
  numeric third = GiNaC::ex_to<numeric>(
      _dF[2].subs(GiNaC::lst{_x == P.x(), _y == P.y(), _z == P.z()}).evalf());

  return Vector(first, second, third);
}

// returns some tangent vector at the point
Vector Function::get_tangent_at_point(const Point &P) const {
  return (get_gradient_at_point(P)).get_any_perpendicular();
}
numeric Function::eval_at_point(const Point &P) const {
  return GiNaC::ex_to<numeric>(
      _F.subs(GiNaC::lst{_x == P.x(), _y == P.y(), _z == P.z()}).evalf());
}
// inside is defined as {P : F(P) < 0}
bool Function::is_inside(const Point &P) const {
  return (eval_at_point(P) <= -kEps);
}
bool Function::is_on(const Point &P) const {
  return abs(eval_at_point(P)) < kEps;
}
bool Function::is_outside(const Point &P) const {
  return (eval_at_point(P) >= kEps);
}

Vector Function::outside_normal(const Triangle &T, const numeric e_size) const {
  // std::cout << "Triangle " << T << endl;
  Vector normal = T.get_normal();
  Point A = numeric(1) / numeric(3) * (T.A() + T.B() + T.C());
  if (is_inside(Point(A, numeric(0.1) * e_size * normal)) &&
      !is_inside(Point(A, -numeric(0.1) * e_size * normal))) {
    // std::cout << "case1" << endl;
    return numeric(-1) * normal;
  } else if (!is_inside(Point(A, numeric(0.1) * e_size * normal)) &&
             is_inside(Point(A, -numeric(0.1) * e_size * normal))) {
    // std::cout << "case2" << endl;
    return normal;
  } else if (is_inside(Point(A, numeric(0.1) * e_size * normal)) &&
             is_inside(Point(A, -numeric(0.1) * e_size * normal))) {
    // std::cout << "case3" << endl;
    return outside_normal(T, e_size * 1.1);
  } else if (!is_inside(Point(A, numeric(0.1) * e_size * normal)) &&
             !is_inside(Point(A, -numeric(0.1) * e_size * normal))) {
    // std::cout << "case4" << endl;
    return outside_normal(T, e_size / 1.1);
  }
  assertm(false, "Should not get here!");
  return normal;
}

numeric Function::substitute(GiNaC::ex il) const {
  return GiNaC::ex_to<numeric>(get_function().subs(il).evalf());
}