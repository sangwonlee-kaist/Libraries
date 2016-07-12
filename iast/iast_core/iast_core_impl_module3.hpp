void
iast_core::setTotalLoading(double totalLoading)
    {
    mTotalLoading = totalLoading;
    }

struct iast_core::Result3
    {
    double totalPressure;
    vec gasFractions;
    vec spreadingPressures;
    vec particularPressures;
    };

std::ostream&
operator << (std::ostream& os, iast_core::Result3& result)
    {
    os << std::setw(15) << result.totalPressure;

    for (const auto& ld : result.gasFractions)
        os << std::setw(15) << ld;

    for (const auto& sp : result.spreadingPressures)
        os << std::setw(15) << sp;

    for (const auto& pp : result.particularPressures)
        os << std::setw(15) << pp;

    return os;
    }

iast_core::Result3
iast_core::calculate3(OBJECTIVE objective)
    {
    // New version, fraction based root finding.
    using RootFinder = root_finder;

    auto  T = m_temperature;
    auto  N = mTotalLoading;
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
            rootFinder.add_eqn([i, j, &f_i, &T](const RootFinder::vec& p)
                {
                return f_i[i](T, p(i)) - f_i[j](T, p(j));
                });
            }
        else if (objective == POW1)
            {
            rootFinder.add_eqn([i, j, &f_i, &T](const RootFinder::vec& p)
                {
                return f_i[i](T, p(i)) / f_i[j](T, p(j)) - 1.0;
                });
            }
        else if (objective == POW2)
            {
            rootFinder.add_eqn([i, j, &f_i, &T](const RootFinder::vec& p)
                {
                return std::pow(f_i[i](T, p(i)) / f_i[j](T, p(j)), 2) - 1.0;
                });
            }
        else if (objective == POW3)
            {
            rootFinder.add_eqn([i, j, &f_i, &T](const RootFinder::vec& p)
                {
                return std::pow(f_i[i](T, p(i)) / f_i[j](T, p(j)), 3) - 1.0;
                });
            }
        else if (objective == POW4)
            {
            rootFinder.add_eqn([i, j, &f_i, &T](const RootFinder::vec& p)
                {
                return std::pow(f_i[i](T, p(i)) / f_i[j](T, p(j)), 4) - 1.0;
                });
            }
        else
            {
            throw std::invalid_argument
                {"iast_core::calculate3() invalid method index"};
            }
        }

    rootFinder.add_eqn([&T, &x_i, &n_i, &N, n_comp](const RootFinder::vec& p)
        {
        RootFinder::vec ni;
        ni.resize(n_comp);
        for (size_t i = 0; i < p.n_elem; ++i)
            ni(i) = n_i[i](T, p(i));

        return N * arma::sum(x_i / ni) - 1.0;
        });

    RootFinder::vec po_i;
    po_i.resize(n_comp);

    // Initial particular pressure guess.
    for (size_t i = 0; i < n_comp; ++i)
        {
        double ni = N * x_i(i);
        double oldP = 1.0; // 1 bar.
        double newP = 0.0;
        for (int iter = 0; iter < 100; ++iter)
            {
            newP = oldP / n_i[i](T, oldP) * ni;
            if (std::abs(1.0 - oldP / newP) < 1.e-3)
                break;
            oldP = newP;
            }
        po_i(i) = newP;
        }

    rootFinder.set_initial_guess(po_i);
    po_i = rootFinder.solve();

    // Get total pressure.
    double totalPressure = arma::dot(po_i, x_i);

    vec ys;
    ys.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        ys[i] = po_i(i) * x_i(i) / totalPressure;

    vec sp;
    sp.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        sp[i] = f_i[i](T, po_i(i));

    vec parp;
    parp.resize(n_comp);

    for (size_t i = 0; i <= last_i; ++i)
        parp[i] = po_i(i);

    return Result3 {totalPressure, ys, sp, parp};
    }
