#include "dsl_isotherm.hpp"

#include <cmath>
#include <sstream>

DslIsotherm::DslIsotherm(double q1, double k1, double q2, double k2) :
    iso1 {q1, k1}, iso2 {q2, k2}, params {q1, k1, q2, k2}
    {

    }

double
DslIsotherm::loading(double P) const
    {
    return iso1.loading(P) + iso2.loading(P);
    }

double
DslIsotherm::spressure(double P) const
    {
    return iso1.spressure(P) + iso2.spressure(P);
    }

std::string
DslIsotherm::getInfoString() const
    {
    std::stringstream ss;

    ss << "[Dual Site Langmuir Isotherm]\n";
    ss << "[Parameters] " <<
          "q1 = " << params[0] << ", " <<
          "k1 = " << params[1] << ", " <<
          "q2 = " << params[2] << ", " <<
          "k2 = " << params[3];

    return ss.str();
    }
