#include "langmuir_isotherm.hpp"

#include <cmath>
#include <sstream>

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
    ss << "langmuir" << endl;
    for (const auto& e : getParameters())
        cout << e.first << "  " << e.second << endl;

    return ss.str();
    }

Isotherm::ParameterType
LangmuirIsotherm::getParameters() const
    {
    ParameterType params;
    params["q"] = mQ;
    params["k"] = mK;

    return params;
    } 
