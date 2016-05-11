#pragma once
#ifndef ISOTHERM_BASE_IMPL_HPP
#define ISOTHERM_BASE_IMPL_HPP

isotherm_base::func_t
isotherm_base::get_loading()
    {
    return [&](real_t T, real_t P) -> real_t
        {
        return loading(T, P);    
        };
    }

isotherm_base::func_t 
isotherm_base::get_spreading_pressure()
    {
    return [&](real_t T, real_t P) -> real_t 
        {
        return spreading_pressure(T, P);
        };
    }

#endif