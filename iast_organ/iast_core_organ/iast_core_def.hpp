#pragma once
#ifndef IAST_CORE_DEF_HPP
#define IAST_CORE_DEF_HPP

class iast_core
    {
public:
    typedef double real_t;
    typedef std::valarray<real_t> vec;
    typedef std::valarray<std::array<real_t, 2>> sp_mat;
    typedef std::vector<isotherm_base::func_t> func_vec;

    void add_isotherm(isotherm_base& iso);

    void set_temperature(real_t temper);    
    void set_pressure(real_t pressure);
    void set_composition(const vec& ys);
    void set_initial_guess(const vec& loading_fractions); 
    int size();
    
    class result;
    
    result calculate();
        
    
    
protected:
    bool     m_has_initial_guess {false};
    real_t   m_temperature;
    real_t   m_pressure;
    vec      m_composition;  
    vec      m_initial_guess;
    func_vec m_loadings;
    func_vec m_spreading_pressures;
    
    };

#endif
