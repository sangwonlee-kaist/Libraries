#pragma once
#ifndef IAST_CORE_IMPL_HPP
#define IAST_CORE_IMPL_HPP

void 
iast_core::add_isotherm(isotherm_base& iso)
    {
    m_loadings.push_back(iso.get_loading());
    m_spreading_pressures.push_back(iso.get_spreading_pressure());
    }

void 
iast_core::set_temperature(real_t temper)
    {
    m_temperature = temper;
    }
    
void
iast_core::set_pressure(real_t pressure)
    {
    m_pressure = pressure;
    }
    
void 
iast_core::set_composition(const vec& ys)
    {
    m_composition = (ys / ys.sum());
    }
    
void 
iast_core::set_initial_guess(const vec& loading_fractions)
    {
    m_initial_guess = loading_fractions /  loading_fractions.sum(); 
    m_has_initial_guess = true;
    }

class iast_core::result
    {
public:        
    result(const vec& in_loadings, const size_t& in_cycle)
        :
        loadings {in_loadings},
        cycle    {in_cycle}
        {
        // nothing.    
        }
    
    vec 
    get_loadings()
        {
        return loadings;      
        }
        
    real_t
    get_total_loading()
        {
        return loadings.sum();
        }

    size_t
    get_cycle()
        {
        return cycle;
        }
        
private:
    vec loadings;
	size_t cycle;	
    };
    
iast_core::result
iast_core::calculate()
    {
    auto  T   = m_temperature;
    auto  P   = m_pressure;
    auto& n_i = m_loadings;
    auto& f_i = m_spreading_pressures;
//    auto& y_i = m_composition;
    auto  n_comp = m_composition.size();
    auto  last_i = n_comp - 1;
    
    root_finder::vec po_i;
    po_i.resize(n_comp);
    
    root_finder::vec y_i;
    y_i.resize(n_comp);

    for (size_t i {}; i < n_comp; ++i)
        {
        y_i(i) = m_composition[i];
        }

    root_finder::vec x_i;
    x_i.resize(n_comp);
    
    if (m_has_initial_guess)
        {
//        x_i = m_initial_guess;
        for (size_t i {}; i < n_comp; ++i)
            {
            x_i(i) = m_initial_guess[i];
            }
        m_has_initial_guess = false;
        }
    else 
        {
        // guess adsorbed phase fraction.
        for (size_t i {}; i < n_comp; ++i)
            {
            x_i[i] = n_i[i](T, P * y_i[i]);
            }
//        x_i /= x_i.sum();            
        x_i = arma::normalise(x_i);        
        }
        
    po_i = P * y_i / x_i;
   
    // find root.
    root_finder rf;
    for (size_t i {0}; i < last_i; ++i)
        {
        rf.add_eqn([i, last_i, &f_i, &T](const root_finder::vec& p)
            {
            return std::sin(f_i[i](T, p(i)) - f_i[last_i](T, p(last_i)));
            });
        }
 
    rf.add_eqn([&P, &y_i](const root_finder::vec& p)
        {
        return P * arma::sum(y_i / p) - 1.0;        
        });

    rf.set_initial_guess(po_i);
    po_i = rf.solve();
    size_t iters = rf.get_iterations();
    x_i = P * y_i / po_i;
        
    double nt {};
            
    for (size_t i {}; i <= last_i; ++i)
        {
        nt += x_i(i) / n_i[i](T, po_i[i]);
        }
    nt = 1.0 / nt;
   
    vec xs;
    xs.resize(n_comp);

    for (size_t i {}; i <= last_i; ++i)
        {
        xs[i] = x_i(i);
        }
 
    return result {nt * xs, iters};
    }

#endif
