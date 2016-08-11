#pragma once

class Isotherm
    {
public:
    virtual ~Isotherm()

    virtual double loading(double T, double P) = 0;
    virtual double spressure(double T, double P) = 0;
    };
