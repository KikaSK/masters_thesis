#pragma once

#include "constants.h"

void set_epsylon(GiNaC::numeric edge_size) {
  kEps = edge_size * GiNaC::numeric(10e-6);
}