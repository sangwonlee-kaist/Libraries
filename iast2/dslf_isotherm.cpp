#include "dslf_isotherm.hpp"

#include <sstream>

DslfIsotherm::DslfIsotherm(double q1, double k1, double n1,
                           double q2, double k2, double n2) :
    iso1 {q1, k1, n1},
    iso2 {q2, k2, n2},
    params {q1, k1, n1, q2, k2, n2}
    {

    }

double
DslfIsotherm::loading(double P) const
    {
    return iso1.loading(P) + iso2.loading(P);
    }

double
DslfIsotherm::spressure(double P) const
    {
    return iso1.spressure(P) + iso2.spressure(P);
    }

std::string
DslfIsotherm::getInfoString() const
    {
    std::stringstream ss;
    ss << "dslf" << '\n';
    for (const auto& e : getParameters())
        ss << e.first << "  " << e.second << '\n';

    return ss.str();
    }

Isotherm::ParameterType
DslfIsotherm::getParameters() const
    {
    ParameterType paramss;

    paramss["q1"] = params[0];
    paramss["k1"] = params[1];
    paramss["n1"] = params[2];
    paramss["q2"] = params[3];
    paramss["k2"] = params[4];
    paramss["n2"] = params[5];

    return paramss;
    }
