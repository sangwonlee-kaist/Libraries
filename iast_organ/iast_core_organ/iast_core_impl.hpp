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
    for (size_t i {}; i < last_i; ++i)
        {
        rf.add_eqn([i, &f_i, &T](const root_finder::vec& p)
            {
            return f_i[i](T, p(i)) / f_i[i + 1](T, p(i + 1)) - 1; 
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

//iast_core::result
//iast_core::calculate()
//    {
//    auto  T   = m_temperature;
//    auto  P   = m_pressure;
//    auto& n_i = m_loadings;
//    auto& f_i = m_spreading_pressures;
//    auto& y_i = m_composition;
//    auto  n_comp = m_composition.size();
//    auto  last_i = n_comp - 1;
//    
//    vec po_i;
//    po_i.resize(n_comp);
//    
//    vec x_i;
//    x_i.resize(n_comp);
//    
//    vec del_i;
//    del_i.resize(n_comp);
//    
//    if (m_has_initial_guess)
//        {
//        x_i = m_initial_guess;
//        m_has_initial_guess = false;
//        }
//    else 
//        {
//        // guess adsorbed phase fraction.
//        for (size_t i {}; i < n_comp; ++i)
//            {
//            x_i[i] = n_i[i](T, P * y_i[i]);
//            }
//        x_i /= x_i.sum();            
//        }
//        
//    po_i = P * y_i / x_i;
//    
//    sp_mat a_ij;
//    a_ij.resize(n_comp);
//    
//    vec b_i;
//    b_i.resize(n_comp);
//    
//    // eq. (20)
//    auto dfdp = [&](size_t i) -> real_t
//        {
//        return n_i[i](T, po_i[i]) / po_i[i];
//        };
//    
//    // eq. (14) and eq. (28)
//    auto g = [&](size_t i) -> real_t
//        {
//        if (i != last_i)
//            {
//            return f_i[i](T, po_i[i]) - f_i[last_i](T, po_i[last_i]);
//            }
//        else    
//            {
//            return P * (y_i / po_i).sum() - 1.0;    
//            }
//        };
//    
//    // eq. (22)
//    auto a_Ni = [&](size_t i) -> real_t
//        {
//        return -P * y_i[i] / po_i[i] / po_i[i];
//        };
//
//	double alpha {0.01};
//    for (int iter {}; true; ++iter)
//        { 
//        double dfdN {dfdp(last_i)};
//        
//        for (size_t i {}; i < last_i; ++i)
//            {
//            a_ij[i][0] =  dfdp(i);
//            a_ij[i][1] = -dfdN;
//            }
//            
//        a_ij[last_i][0] = 0.0;
//        a_ij[last_i][1] = a_Ni(last_i);
//        
//        for (size_t i {}; i <= last_i; ++i)
//            {
//            b_i[i] = -g(i);
//            }
//        
//        for (size_t i {}; i < last_i; ++i)
//            {
//            a_ij[last_i][1] -= a_Ni(i) / a_ij[i][0] * a_ij[i][1];
//            b_i[last_i]     -= a_Ni(i) / a_ij[i][0] * b_i[i];
//            }
//            
//        del_i[last_i] = b_i[last_i] / a_ij[last_i][1];
//        double xn {del_i[last_i]};
//        
//        for (size_t i {}; i < last_i; ++i)
//            {
//            del_i[i] = b_i[i] - a_ij[i][1] * xn / a_ij[i][0];
//            }
//        
//        del_i *= alpha;
//        po_i  += del_i;
//        
//        if (std::any_of(std::begin(po_i), std::end(po_i), [](double x) {return x < 0;}))
//            {
////            po_i  -= del_i;
////            del_i *= 0.99;
////            po_i  += del_i;
//            auto old_po_i = po_i;              
//            po_i  = std::abs(po_i);
//            del_i = po_i - old_po_i; 
//            }
//     
//        if (not std::all_of(std::begin(po_i), std::end(po_i), 
//            [](double x) {return std::isfinite(x);}))
//            {
//            throw std::runtime_error {"[error] iast is not converged."};
//            } 
//        
//        auto reduce_i = po_i.apply([](double x) {return (x > 1.0 ? x : 1.0);});  
//        double error = std::abs(del_i / reduce_i).sum();
////        double error = std::abs(del_i / po_i).sum();   
//
////        if (iter % 10000 == 0)
//          if (false)
//            {
//            for (auto& po : po_i)
//                {
//                std::cout << std::setw(15) << po;
//                }
//            for (auto& del : del_i)
//                {
//                std::cout << std::setw(15) << del;
//                }
//            std::cout << std::setw(15) << error << std::endl;
//            }
//
//        if (error < 1.e-10)
//            {
//            vec func_test;
//            func_test.resize(last_i + 1);
//
//            for (size_t i {}; i <= last_i; ++i)
//                {
//                func_test[i] = std::abs(g(i)); 
//                }
//
//            if (func_test.min() > 1.e-7)
//                {
//                for (auto fn : func_test)
//                    {
//                    std::cout << std::setw(15) << fn;
//                    }
//                std::cout << std::endl;
//                throw std::runtime_error {"[error] objective vector is not zero."};
//                }
//        
//            x_i = P * y_i / po_i;
//            
//            double nt {};
//            
//            for (size_t i {}; i <= last_i; ++i)
//                {
//                nt += x_i[i] / n_i[i](T, po_i[i]);
//                }
//            nt = 1.0 / nt;
//            
//            return result {nt * x_i, static_cast<size_t>(iter)};
//            //break;
//            }
//        
//        if (iter > 1e6)
//            {
//            throw std::runtime_error {"[error] iast exceed maximum iteration"};     
//            }
//        } // iter
//    }
//
#endif
