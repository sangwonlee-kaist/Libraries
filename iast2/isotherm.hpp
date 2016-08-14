#pragma once

class Isotherm
    {
public:
    virtual ~Isotherm() = default;

    virtual double loading(double P) = 0;
    virtual double spressure(double P) = 0;
    };
