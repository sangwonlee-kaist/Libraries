#pragma once
#ifndef ISOTHERM_BASE_DEF_HPP
#define ISOTHERM_BASE_DEF_HPP

class isotherm_base
    {
public:
    typedef double real_t;
    typedef std::function<real_t(real_t, real_t)> func_t;

    func_t get_loading() const;
    func_t get_spreading_pressure() const;
    
    // units: temperature = K and pressure = bar. 
    virtual real_t loading(real_t temper, real_t press) const          = 0;
    virtual real_t spreading_pressure(real_t temper, real_t press) const = 0;
    };

#endif
