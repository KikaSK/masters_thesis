#pragma once

#include <cassert>
#include <iostream>

#define assertm(exp, msg) assert(((void)msg, exp))
