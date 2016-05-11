#pragma once

// check c++11.
#if (__cplusplus >= 201103L)
    #undef  SIMPLEX_USE_CXX11
    #define SIMPLEX_USE_CXX11
#endif

#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdexcept>

#ifdef SIMPLEX_USE_CXX11
    #include <random>
    #include <chrono>
#else
    #include <ctime>
#endif

#include "simplex/simplex_def.hpp"
#include "simplex/simplex_impl.hpp"
