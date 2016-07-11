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
iast_core::set_composition(const std::vector<double>& ys)
    {
    vec temp (ys.data(), ys.size());
    this->set_composition(temp);
    }

void
iast_core::set_initial_guess(const vec& particular_pressures)
    {
    m_initial_guess = particular_pressures;
    m_has_initial_guess = true;
    }

int
iast_core::size()
    {
    return static_cast<int>(m_loadings.size());
    }

class iast_core::result
    {
public:
    result(const vec& in_loadings,
           const vec& in_spreading_pressures,
           const vec& in_particular_pressure,
           const size_t& in_cycle)
        :
        loadings             {in_loadings},
        spreading_pressures  {in_spreading_pressures},
        particular_pressures {in_particular_pressure},
        cycle                {in_cycle}
        {
        // nothing.
        }

    vec
    get_loadings()
        {
        return loadings;
        }

    vec
    get_spreading_pressures_for_checking()
        {
        return spreading_pressures;
        }

    vec
    get_particular_pressures()
        {
        return particular_pressures;
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
    vec spreading_pressures;
    vec particular_pressures;
    size_t cycle;
    };

std::ostream&
operator << (std::ostream& os, iast_core::result& result)
    {
    for (const auto& ld : result.get_loadings())
        os << std::setw(15) << ld;

    for (const auto& sp : result.get_spreading_pressures_for_checking())
        os << std::setw(15) << sp;

    for (const auto& pp : result.get_particular_pressures())
        os << std::setw(15) << pp;

    return os;
    }

iast_core::result
iast_core::calculate(OBJECTIVE objective)
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
        // Put initial guess.
        // po_i = m_initial_guess;
        for (size_t i = 0; i < n_comp; ++i)
            {
            po_i(i) = m_initial_guess[i];
            }
        x_i = P * y_i / po_i;
        // Turn off initial guess flag.
        m_has_initial_guess = false;
        }
    else
        {
        // guess adsorbed phase fraction.
        for (size_t i {}; i < n_comp; ++i)
            {
            x_i[i] = n_i[i](T, P * y_i[i]);
            }
        x_i = arma::normalise(x_i);
        po_i = P * y_i / x_i;
        }

    // find root.
    root_finder rf;
    for (size_t i = 0; i <= last_i; ++i)
        {
        size_t pivotIndex = 0;
        if (i == pivotIndex)
            continue;

        if (objective == DIFF)
            {
            rf.add_eqn([i, last_i, &f_i, &T, pivotIndex](const root_finder::vec& p)
                {
                return f_i[i](T, p(i)) - f_i[pivotIndex](T, p(pivotIndex));
                });
            }
        else if (objective == POW1)
            {
            rf.add_eqn([i, last_i, &f_i, &T, pivotIndex](const root_finder::vec& p)
                {
                return  f_i[i](T, p(i)) / f_i[pivotIndex](T, p(pivotIndex)) - 1.0;
                });
            }
        else if (objective == POW2)
            {
            rf.add_eqn([i, last_i, &f_i, &T, pivotIndex](const root_finder::vec& p)
                {
                return std::pow(f_i[i](T, p(i)) / f_i[pivotIndex](T, p(pivotIndex)), 2) - 1.0;
                });
            }
        else if (objective == POW3)
            {
            rf.add_eqn([i, last_i, &f_i, &T, pivotIndex](const root_finder::vec& p)
                {
                return std::pow(f_i[i](T, p(i)) / f_i[pivotIndex](T, p(pivotIndex)), 3) - 1.0;
                });
            }
        else if (objective == POW4)
            {
            rf.add_eqn([i, last_i, &f_i, &T, pivotIndex](const root_finder::vec& p)
                {
                return std::pow(f_i[i](T, p(i)) / f_i[pivotIndex](T, p(pivotIndex)), 4) - 1.0;
                });
            }
        }

    rf.add_eqn([&P, &y_i](const root_finder::vec& p)
        {
        return P * arma::sum(y_i / p) - 1.0;
        //return std::log(P * arma::sum(y_i / p));
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

    vec sp;
    sp.resize(n_comp);

    for (size_t i {}; i <= last_i; ++i)
        {
        sp[i] = f_i[i](T, po_i(i));
        }

    vec parp;
    parp.resize(n_comp);

    for (size_t i {}; i <= last_i; ++i)
        {
        parp[i] = po_i(i);
        }

    return result {nt * xs, sp, parp, iters};
    }

iast_core::result
iast_core::calculate1(OBJECTIVE objective)
    {
    // New version, fraction based root finding.
    using RootFinder = root_finder;

    auto  T   = m_temperature;
    auto  P   = m_pressure;
    auto& n_i = m_loadings;
    auto& f_i = m_spreading_pressures;
    auto  n_comp = m_composition.size();
    auto  last_i = n_comp - 1;

    RootFinder::vec y_i;
    y_i.resize(n_comp);

    for (size_t i = 0; i < n_comp; ++i)
        y_i(i) = m_composition[i];

    // finding root algorithm.
    RootFinder rootFinder;

    // Set objective function.
    for (size_t i = 0; i <= last_i; ++i)
        {
        size_t pivotIndex = 0;
        size_t j = pivotIndex;
        if (i == pivotIndex)
            continue;

        if (objective == DIFF)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &y_i](const RootFinder::vec& x)
                {
                return f_i[i](T, P * y_i(i) / x(i)) -
                       f_i[j](T, P * y_i(j) / x(j));
                });
            }
        else if (objective == POW1)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &y_i](const RootFinder::vec& x)
                {
                return f_i[i](T, P * y_i(i) / x(i)) /
                       f_i[j](T, P * y_i(j) / x(j)) - 1.0;
                });
            }
        else if (objective == POW2)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &y_i](const RootFinder::vec& x)
                {
                return std::pow(f_i[i](T, P * y_i(i) / x(i)) /
                                f_i[j](T, P * y_i(j) / x(j)), 2) - 1.0;
                });
            }
        else if (objective == POW3)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &y_i](const RootFinder::vec& x)
                {
                return std::pow(f_i[i](T, P * y_i(i) / x(i)) /
                                f_i[j](T, P * y_i(j) / x(j)), 3) - 1.0;
                });
            }
        else if (objective == POW4)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &y_i](const RootFinder::vec& x)
                {
                return std::pow(f_i[i](T, P * y_i(i) / x(i)) /
                                f_i[j](T, P * y_i(j) / x(j)), 4) - 1.0;
                });
            }
        else
            {
            throw std::invalid_argument
                {"iast_core::calculate1() invalid method index"};
            }
        }

    rootFinder.add_eqn([](const RootFinder::vec& x)
        {
        return arma::sum(x) - 1.0;
        });

    RootFinder::vec x_i;
    x_i.resize(n_comp);

    // guess adsorbed phase fraction.
    for (size_t i = 0; i < n_comp; ++i)
        x_i[i] = n_i[i](T, P * y_i[i]);
    x_i = arma::normalise(x_i);

    rootFinder.set_initial_guess(x_i);
    x_i = rootFinder.solve();

    size_t iters = rootFinder.get_iterations();
    // Get particular pressures.
    RootFinder::vec po_i;
    po_i.resize(n_comp);
    po_i = P * y_i / x_i;

    double nt = 0.0;
    for (size_t i = 0; i <= last_i; ++i)
        nt += x_i(i) / n_i[i](T, po_i[i]);
    nt = 1.0 / nt;

    vec xs;
    xs.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        xs[i] = x_i(i);

    vec sp;
    sp.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        sp[i] = f_i[i](T, po_i(i));

    vec parp;
    parp.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        parp[i] = po_i(i);

    return result {nt * xs, sp, parp, iters};
    }

#endif
