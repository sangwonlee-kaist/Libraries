#pragma once

#include <cmath>
#include <sstream>
#include "isotherm.hpp"
#include "lf_isotherm.hpp"

class DslfIsotherm : public Isotherm
    {
public:
    DslfIsotherm(double q1, double k1, double n1,
                 double q2, double k2, double n2);

    virtual ~DslfIsotherm() = default;

    virtual double loading(double P) const override;
    virtual double spressure(double P) const override;

    virtual std::string getInfoString() const override;
private:
    LfIsotherm iso1;
    LfIsotherm iso2;
    double params[6];
    };

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

    ss << "[Dual Site Langmuir Freundlich Isotherm]\n";
    ss << "[Parameters] " <<
          "q1 = " << params[0] << ", " <<
          "k1 = " << params[1] << ", " <<
          "n1 = " << params[2] << ", " <<
          "q2 = " << params[3] << ", " <<
          "k2 = " << params[4] << ", " <<
          "n2 = " << params[5];

    return ss.str();
    }
