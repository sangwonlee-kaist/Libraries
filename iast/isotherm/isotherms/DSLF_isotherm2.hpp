#pragma once

#include <cmath>

class DSLF_isotherm2 : public isotherm_base
    {
public:
    DSLF_isotherm2(real_t in_qm1, real_t in_qm2, real_t in_b1, real_t in_b2, real_t in_n1, real_t in_n2)
        :
        qm1 {in_qm1},
        qm2 {in_qm2},
        b1  {in_b1},
        b2  {in_b2},
        n1  {in_n1},
        n2  {in_n2}
        {
        // do nothing.
        }

    real_t
    loading(real_t T, real_t P)
        override
        {
        double p = P;
        if (P ==  0.0)
            {
            return 0.0;
            }

        if (P < 0.0)
            {
//            return loading(T, -P);
            return 0.0;
            }

        double p_n1 {std::pow(p, n1)};
        double p_n2 {std::pow(p, n2)};

        return qm1 * b1 * p_n1 / (1.0 + b1 * p_n1) + qm2 * b2 * p_n2 / (1.0 + b2 * p_n2);
        }

    real_t
    spreading_pressure(real_t T, real_t P)
        override
        {
        double p = P;
        if (P ==  0.0)
            {
            return 0.0;
            }

        if (P < 0.0)
            {
//            return spreading_pressure(T, -P);
            return 0.0;
            }

        double p_n1 {std::pow(p, n1)};
        double p_n2 {std::pow(p, n2)};

        return qm1 / n1 * std::log(1.0 + b1 * p_n1) + qm2 / n2 * std::log(1.0 + b2 * p_n2);
        }

private:
    double qm1;
    double qm2;
    double b1;
    double b2;
    double n1;
    double n2;

    };
