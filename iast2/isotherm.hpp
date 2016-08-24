#pragma once

#include <string>

class Isotherm
    {
public:
    virtual ~Isotherm() = default;

    virtual double loading(double P) const = 0;
    virtual double spressure(double P) const = 0;

    virtual std::string getInfoString() const = 0;
    };
