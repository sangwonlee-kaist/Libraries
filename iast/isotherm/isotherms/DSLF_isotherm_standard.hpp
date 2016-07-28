#pragma once

class DslfIsotherm : public isotherm_base
    {
public:
    DslfIsotherm(double q1_, double k1_, double n1_,
                 double q2_, double k2_, double n2_) :
        q1 (q1_), k1 (k1_), n1 (n1_),
        q2 (q2_), k2 (k2_), n2 (n2_)
        {
        // do nothing.
        }

    real_t
    loading(real_t T, real_t P) const override
        {
        if (P <= 0.0)
            return 0.0;

        double term1 = std::pow(k1 * P, n1);
        double term2 = std::pow(k2 * P, n2);

        return q1 * term1 / (1.0 + term1) + q2 * term2 / (1.0 + term2);
        }

    real_t
    spreading_pressure(real_t T, real_t P) const override
        {
        if (P <=  0.0)
            return 0.0;

        double term1 = std::pow(k1 * P, n1);
        double term2 = std::pow(k2 * P, n2);

        return q1 * std::log(1.0 + term1) / n1 + q2 * std::log(1.0 + term2) / n2;
        }
private:
    double q1;
    double q2;
    double k1;
    double k2;
    double n1;
    double n2;
    };
