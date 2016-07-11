#pragma once
#ifndef IAST_CORE_DEF_HPP
#define IAST_CORE_DEF_HPP

enum OBJECTIVE
    {
    DIFF = 0,
    POW1,
    POW2,
    POW3,
    POW4
    };

class iast_core
    {
public:
    typedef double real_t;
    typedef std::valarray<real_t> vec;
    typedef std::vector<isotherm_base::func_t> func_vec;

    void add_isotherm(isotherm_base& iso);

    void set_temperature(real_t temper);
    void set_pressure(real_t pressure);
    void set_composition(const vec& ys);
    void set_composition(const std::vector<double>& ys);
    void set_initial_guess(const vec& particular_pressures);
    int size();

    class result;

    result calculate(OBJECTIVE objective = DIFF);
protected:
    bool     m_has_initial_guess {false};
    real_t   m_temperature;
    real_t   m_pressure;
    vec      m_composition;
    vec      m_initial_guess;
    func_vec m_loadings;
    func_vec m_spreading_pressures;

public:
    // For fixed Temperature...
    // 1) Find x_i at given P and y_i
    result calculate1(OBJECTIVE objective = DIFF);
    };

#endif
