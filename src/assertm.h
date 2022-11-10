#ifndef ASSERTM_H
#define ASSERTM_H

#include <cassert>
#include <iostream>

#define assertm(exp, msg) assert(((void)msg, exp))

#endif