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
    ss << "lf" << endl;
    for (const auto& e : getParameters())
        cout << e.first << "  " << e.second << endl;

    return ss.str();
    }


Isotherm::ParameterType
LfIsotherm::getParameters() const
    {
    ParameterType params;

    params["q"] = mQ;
    params["k"] = mK;
    params["n"] = mN;

    return params;
    }
