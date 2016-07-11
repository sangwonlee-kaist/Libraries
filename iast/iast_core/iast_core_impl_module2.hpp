struct iast_core::Result2
    {
    double totalUptake;
    vec gasFractions;
    vec spreadingPressures;
    vec particularPressures;
    };

std::ostream&
operator << (std::ostream& os, iast_core::Result2& result)
    {
    os << std::setw(15) << result.totalUptake;

    for (const auto& ld : result.gasFractions)
        os << std::setw(15) << ld;

    for (const auto& sp : result.spreadingPressures)
        os << std::setw(15) << sp;

    for (const auto& pp : result.particularPressures)
        os << std::setw(15) << pp;

    return os;
    }

iast_core::Result2
iast_core::calculate2(OBJECTIVE objective)
    {
    // New version, fraction based root finding.
    using RootFinder = root_finder;

    auto  T   = m_temperature;
    auto  P   = m_pressure;
    auto& n_i = m_loadings;
    auto& f_i = m_spreading_pressures;
    auto  n_comp = m_composition.size();
    auto  last_i = n_comp - 1;

    RootFinder::vec x_i;
    x_i.resize(n_comp);

    for (size_t i = 0; i < n_comp; ++i)
        x_i(i) = m_composition[i];

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
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &x_i](const RootFinder::vec& y)
                {
                return f_i[i](T, P * y(i) / x_i(i)) -
                       f_i[j](T, P * y(j) / x_i(j));
                });
            }
        else if (objective == POW1)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &x_i](const RootFinder::vec& y)
                {
                return f_i[i](T, P * y(i) / x_i(i)) /
                       f_i[j](T, P * y(j) / x_i(j)) - 1.0;
                });
            }
        else if (objective == POW2)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &x_i](const RootFinder::vec& y)
                {
                return std::pow(f_i[i](T, P * y(i) / x_i(i)) /
                                f_i[j](T, P * y(j) / x_i(j)), 2) - 1.0;
                });
            }
        else if (objective == POW3)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &x_i](const RootFinder::vec& y)
                {
                return std::pow(f_i[i](T, P * y(i) / x_i(i)) /
                                f_i[j](T, P * y(j) / x_i(j)), 3) - 1.0;
                });
            }
        else if (objective == POW4)
            {
            rootFinder.add_eqn([i, j, &f_i, &T, &P, &x_i](const RootFinder::vec& y)
                {
                return std::pow(f_i[i](T, P * y(i) / x_i(i)) /
                                f_i[j](T, P * y(j) / x_i(j)), 4) - 1.0;
                });
            }
        else
            {
            throw std::invalid_argument
                {"iast_core::calculate2() invalid method index"};
            }
        }

    rootFinder.add_eqn([](const RootFinder::vec& y)
        {
        return arma::sum(y) - 1.0;
        });

    RootFinder::vec y_i;
    y_i.resize(n_comp);

    // guess gas phase fraction.
    for (size_t i = 0; i < n_comp; ++i)
        y_i[i] = x_i[i];

    rootFinder.set_initial_guess(y_i);
    y_i = rootFinder.solve();

    size_t iters = rootFinder.get_iterations();
    // Get particular pressures.
    RootFinder::vec po_i;
    po_i.resize(n_comp);
    po_i = P * y_i / x_i;

    double nt = 0.0;
    for (size_t i = 0; i <= last_i; ++i)
        nt += x_i(i) / n_i[i](T, po_i[i]);
    nt = 1.0 / nt;

    vec ys;
    ys.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        ys[i] = y_i(i);

    vec sp;
    sp.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        sp[i] = f_i[i](T, po_i(i));

    vec parp;
    parp.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        parp[i] = po_i(i);

    return Result2 {nt, ys, sp, parp};
    }
