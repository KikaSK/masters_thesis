#include "triangle.h"

#include <iostream>

Triangle::Triangle(Point A, Point B, Point C)
    : _A(A), _B(B), _C(C), _AB(A, B), _BC(B, C), _CA(C, A) {
  assertm(A != B && B != C && A != C,
          "Trying to make triangle from the same points!");

  // assertm(is_triangle(), "Creating non-valid triangle!");
  // if (!(_e1.get_length() + _e2.get_length() - _e3.get_length() >
  // -0.001  &&
  //_e2.get_length() + _e3.get_length() - _e1.get_length() > -0.001 &&
  //_e1.get_length() + _e3.get_length() - _e2.get_length() > -0.001 ))

  // throw "Given points don't represent a Triangle";
};

Point Triangle::A() const { return _A; }
Point Triangle::B() const { return _B; }
Point Triangle::C() const { return _C; }

Edge Triangle::AB() const { return _AB; }
Edge Triangle::BC() const { return _BC; }
Edge Triangle::CA() const { return _CA; }

Point Triangle::get_gravity_center() const {
  assertm(is_triangle(), "Getting gravity center of non valid triangle!");

  auto x_coord = (_A.x() + _B.x() + _C.x()) / 3;
  auto y_coord = (_A.y() + _B.y() + _C.y()) / 3;
  auto z_coord = (_A.z() + _B.z() + _C.z()) / 3;
  return Point(x_coord, y_coord, z_coord);
}
bool Triangle::is_triangle() const {
  // testing triangle inequality
  auto a = _AB.get_length();
  auto b = _BC.get_length();
  auto c = _CA.get_length();
  return (a + b > c + kEps * kEps) && (b + c > a + kEps * kEps) &&
         (a + c > b + kEps * kEps);
}

// returns unit normal vector of triangle
Vector Triangle::get_normal() const {
  assertm(is_triangle(), "Getting normal of non valid triangle!");
  Vector AB(_A, _B);
  Vector AC(_A, _C);

  assertm(!(AB ^ AC).is_zero(), "Normal is zero vector!");
  return (AB ^ AC).unit();
}
// using formula from webpage:
// https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
Point Triangle::get_circumcenter() const {
  assertm(is_triangle(), "Getting circumcenter of non valid triangle!");

  Vector AB(_A, _B);
  Vector AC(_A, _C);

  assertm(_B != _C, "Same edges in triangle.");
  assertm(((AB ^ AC).get_length()) > kEps * kEps,
          "Edges of triangle lineary dependent!");
  Point circumcenter(_A, (AC.get_length_squared() * ((AB ^ AC) ^ AB) +
                          AB.get_length_squared() * ((AC ^ AB) ^ AC)) /
                             (2 * (AB ^ AC).get_length_squared()));
  numeric l1 = Vector(_A, circumcenter).get_length();
  numeric l2 = Vector(_B, circumcenter).get_length();
  numeric l3 = Vector(_C, circumcenter).get_length();
  assertm(abs(l1 - l2) < kEps && abs(l1 - l3) < kEps && abs(l2 - l3) < kEps,
          "Wrong circumcenter!");
  return circumcenter;
}

// finds out whether point P is inside triangle T
bool Triangle::is_in_triangle(Point P) const {
  // using barycentric algorithm from:
  // https://blackpawn.com/texts/pointinpoly/default.html

  // std::cout << endl << endl << "In the is_in_triangle algorithm!" << endl;
  // std::cout<< "Points: " << endl << "A: " << _A << " B: " <<_B << " C: " <<
  // _C << " P: " << P <<endl;
  // std::cout << "triangle " << *this << endl;
  assertm(is_triangle(), "Non valid triangle!");

  // vectors
  Vector v0(_A, _C);
  Vector v1(_A, _B);
  Vector v2(_A, P);

  // std::cout<< "Vectors: " <<endl << "AC: " << v0 << " AB: " << v1 << " AP: "
  // << v2 << endl;

  // dot products
  numeric dot00 = v0 * v0;
  numeric dot01 = v0 * v1;
  numeric dot02 = v0 * v2;
  numeric dot11 = v1 * v1;
  numeric dot12 = v1 * v2;

  numeric denominator = dot00 * dot11 - dot01 * dot01;

  assertm(abs(denominator) > (10e-4) * kEps, "Denominator is zero!");

  numeric u = (dot11 * dot02 - dot01 * dot12) / denominator;
  numeric v = (dot00 * dot12 - dot01 * dot02) / denominator;

  // the point is not in the plane
  if (!(v2 - Vector(u * v0 + v * v1)).is_zero()) return false;

  assertm((v2 - Vector(u * v0 + v * v1)).is_zero(),
          "Not correct barycentric coordinates!");

  // std::cout << endl << "Barycentric coordinates: " << endl << "u: "<<u<<" v:
  // " << v << endl; std::cout << "Leaving in_the_triangle function!" << endl;

  return (u >= 0 && v >= 0 && u + v <= 1);
}