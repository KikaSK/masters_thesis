#include <ginac/ginac.h>

#include <iostream>

#include "mesh_point.h"
#include "point.h"
#include "vector.h"

int main() {
  GiNaC::Digits = 15;

  for (int i = 0; i < 10; ++i) {
    std::cout << i << std::endl;
  }
}
