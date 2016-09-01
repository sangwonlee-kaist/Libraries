#include "henry_isotherm.hpp"

#include <cmath>
#include <sstream>

HenryIsotherm::HenryIsotherm(double k) :
    mK {k}
    {

    }

double
HenryIsotherm::loading(double P) const
    {
    return  mK * P;
    }

double
HenryIsotherm::spressure(double P) const
    {
    return mK * P;
    }

std::string
HenryIsotherm::getInfoString() const
    {
    std::stringstream ss;

    ss << "[Henry Isotherm]\n";
    ss << "[Parameters] k = " << mK;

    return ss.str();
    }
