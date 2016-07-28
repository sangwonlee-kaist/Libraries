#pragma once

class DSL_isotherm : public isotherm_base
    {
public:
    DSL_isotherm(real_t in_q_sat1, real_t in_q_sat2,
                      real_t in_K1,     real_t in_K2)
        :
        q_sat1 {in_q_sat1},
        q_sat2 {in_q_sat2},
        K1  {in_K1},
        K2  {in_K2}
        {
        // do nothing.    
        }
        
    real_t 
    loading(real_t T, real_t P) const override
        {
        return q_sat1 * K1 * P / (1.0 + K1 * P) + q_sat2 * K2 * P / (1.0 + K2 * P);
        }
        
    real_t 
    spreading_pressure(real_t T, real_t P) const override
        {
        return q_sat1 * std::log(1.0 + K1 * P) + q_sat2 * std::log(1.0 + K2 * P);
        }       
private:
    double q_sat1;
    double q_sat2;
    double K1;
    double K2;
    };
