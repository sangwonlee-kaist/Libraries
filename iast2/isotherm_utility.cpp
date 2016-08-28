#include "isotherm.hpp"
#include "isotherm_exception.hpp"

#include "solver.hpp"
#include "solver_factory.hpp"
#include "solver_exception.hpp"

double inverseIsotherm(Isotherm& isotherm, double n)
    {
    auto solver = SolverFactory {}.create("bisection");

    solver->setFunctions({
        [&isotherm, n](const Solver::PointType& p)
            {
            return isotherm.loading(p[0]) - n;
            }
        });

    double p0 = 0.0;
    double p1 = 1.0;

    while (isotherm.loading(p1) < n)
        {
        p1 *= 2.0;

        if (p1 > 10000.0)
            {
            const char* msg {"Given uptake beyonds saturation loading."};
            throw IsothermException {__FILE__, __LINE__, msg};
            }
        }

    try {
        solver->setInitialPoint({p0, p1}).solve();
        }
    catch (SolverException& e)
        {
        const char* msg {"Solver error occurs."};
        throw IsothermException {__FILE__, __LINE__, msg};
        }

    auto pressure = solver->getRootPoint()[0];

    return pressure;
    }
