struct iast_core::Result1
    {
    vec loadings;
    vec spreadingPressures;
    vec particularPressures;
    };

std::ostream&
operator << (std::ostream& os, iast_core::Result1& result)
    {
    for (const auto& ld : result.loadings)
        os << std::setw(15) << ld;

    for (const auto& sp : result.spreadingPressures)
        os << std::setw(15) << sp;

    for (const auto& pp : result.particularPressures)
        os << std::setw(15) << pp;

    return os;
    }

iast_core::Result1
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

    return Result1 {nt * xs, sp, parp};
    }
