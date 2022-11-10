#include <ginac/ginac.h>

#include <iostream>

#include "mesh_point.h"
#include "point.h"
#include "vector.h"

int main() {
  GiNaC::Digits = 15;

  Point P0(0, 0, 0);
  Point Px(1, 0, 0);
  Point Py(0, 1, 0);
  Point Pz(0, 0, 1);

  for (int i = 0; i < 10; ++i) {
    std::cout << i << std::endl;
  }
}
