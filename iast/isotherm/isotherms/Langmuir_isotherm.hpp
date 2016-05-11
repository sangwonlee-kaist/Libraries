#pragma once

#include <cmath>
    
class Langmuir_isotherm : public isotherm_base
    {
public:
    Langmuir_isotherm(real_t in_q_sat, real_t in_K)
        :
        q_sat {in_q_sat},
        K  {in_K}
        {
        // do nothing.    
        }
        
    real_t 
    loading(real_t T, real_t P) 
        override
        {
        return q_sat * K * P / (1.0 + K * P);
        }
        
    real_t 
    spreading_pressure(real_t T, real_t P) 
        override
        {
        return q_sat * std::log(1.0 + K * P);
        }       
        
private:
    double q_sat;
    double K;
    
    };
