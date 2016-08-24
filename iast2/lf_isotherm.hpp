#pragma once

#include <cmath>
#include <sstream>
#include "isotherm.hpp"

class LfIsotherm : public Isotherm
    {
public:
    LfIsotherm(double q, double k, double n);
    virtual ~LfIsotherm() = default;

    virtual double loading(double P) const override;
    virtual double spressure(double P) const override;

    virtual std::string getInfoString() const override;
private:
    double mQ;
    double mK;
    double mN;
    };

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
