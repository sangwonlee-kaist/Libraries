#include "lf_isotherm.hpp"

#include <cmath>
#include <sstream>

LfIsotherm::LfIsotherm(double q, double k, double n) :
    mQ {q}, mK {k}, mN {n}
    {

    }

double
LfIsotherm::loading(double P) const
    {
    double kpn = std::pow(mK * P, mN);
    return mQ * kpn / (1.0 + kpn);
    }

double
LfIsotherm::spressure(double P) const
    {
    double kpn = std::pow(mK * P, mN);
    return mQ / mN * std::log(1.0 + kpn);
    }

std::string
LfIsotherm::getInfoString() const
    {
    std::stringstream ss;

    ss << "[Langmuir Freundlich Isotherm]\n";
    ss << "[Parameters] q = " << mQ << ", k = " << mK << ", n = " << mN;

    return ss.str();
    }
