#pragma once

#include <cmath>

class LF_isotherm : public isotherm_base
    {
public:
    LF_isotherm(real_t in_q_sat, real_t in_K, real_t in_v)
        :
        q_sat {in_q_sat},
        K  {in_K},
        v  {in_v}
        {
        // do nothing.
        }

    real_t
    loading(real_t T, real_t P)
        override
        {
        return q_sat * K * std::pow(P, v) / (1.0 + K * std::pow(P, v));
        }

    real_t
    spreading_pressure(real_t T, real_t P)
        override
        {
        return q_sat / v * std::log(1.0 + K * std::pow(P, v));
        }
private:
    double q_sat;
    double K;
    double v;
    };
