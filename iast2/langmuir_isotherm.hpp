#pragma once

#include <cmath>
#include <sstream>
#include "isotherm.hpp"

class LangmuirIsotherm : public Isotherm
    {
public:
    LangmuirIsotherm(double q, double k);
    virtual ~LangmuirIsotherm() = default;

    virtual double loading(double P) const override;
    virtual double spressure(double P) const override;

    virtual std::string getInfoString() const override;
private:
    double mQ;
    double mK;
    };

LangmuirIsotherm::LangmuirIsotherm(double q, double k) :
    mQ {q}, mK {k}
    {

    }

double
LangmuirIsotherm::loading(double P) const
    {
    return mQ * mK * P / (1.0 + mK * P);
    }

double
LangmuirIsotherm::spressure(double P) const
    {
    return mQ * std::log(1.0 + mK * P);
    }

std::string
LangmuirIsotherm::getInfoString() const
    {
    std::stringstream ss;

    ss << "[Langmuir Isotherm]\n";
    ss << "[Parameters] q = " << mQ << ", k = " << mK;

    return ss.str();
    }
