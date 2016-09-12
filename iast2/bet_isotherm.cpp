#include "bet_isotherm.hpp"

#include <sstream>
#include <cmath>

BetIsotherm::BetIsotherm(double _q, double _k1, double _k2) :
    q {_q}, k1 {_k1}, k2 {_k2}
    {
    }

BetIsotherm::~BetIsotherm()
    {
    }

double
BetIsotherm::loading(double p) const
    {
    return q * k1 * p / (1.0 - k2 * p) / (1.0 + (k1 - k2) * p);
    }

double
BetIsotherm::spressure(double p) const
    {
    return q * std::log((1.0 + (k1 - k2) * p) / (1.0 - k2 * p));
    }

std::string
BetIsotherm::getInfoString() const
    {
    std::stringstream ss;
    ss << "bet" << '\n';
    for (const auto& e : getParameters())
        ss << e.first << "  " << e.second << '\n';

    return ss.str();
    }

Isotherm::ParameterType
BetIsotherm::getParameters() const
    {
    ParameterType params;
    params["q"] = q;
    params["k1"] = k1;
    params["k2"] = k2;

    return params;
    }
