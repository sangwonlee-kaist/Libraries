#include "isotherm_utility.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ciso646>

#include "isotherm_factory.hpp"
#include "isotherm_exception.hpp"

#include "solver.hpp"
#include "solver_factory.hpp"
#include "solver_exception.hpp"

#include "minimizer.hpp"
#include "minimizer_exception.hpp"
#include "simplex.hpp" // ? Factory is needed.

void
readTwoColumns(const std::string& filename,
                     std::vector<double>& x,
                     std::vector<double>& y)
    {
    std::ifstream ifs {filename};

    if (not ifs)
        throw IsothermException {__FILE__, __LINE__, "Invalid filename"};

    readTwoColumns(ifs, x, y);
    }

void
readTwoColumns(std::istream& is,
               std::vector<double>& x,
               std::vector<double>& y)
    {
    x.clear();
    y.clear();

    while (is)
        {
        double xx;
        double yy;

        if (is >> xx >> yy)
            {
            x.push_back(xx);
            y.push_back(yy);
            }
        }
    }

double
inverseIsotherm(Isotherm& isotherm, double n)
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

    double nOld = isotherm.loading(p1);
    while (isotherm.loading(p1) < n)
        {
        p1 *= 2.0;

        double nNew = isotherm.loading(p1);
        //if (p1 > 16384.000001) // 2^14
        if (std::abs(1.0 - nOld / nNew) < 1.0e-4)
            {
            const char* msg {"Given uptake beyonds saturation loading."};
            throw IsothermException {__FILE__, __LINE__, msg};
            }

        nOld = nNew;
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

// IsothermModeler

IsothermModeler::IsothermModeler() : mError {0.0}
    {
    // name : num params
    mIsothermMap["henry"]     = 1;
    mIsothermMap["langmuir"]  = 2;
    mIsothermMap["lf"]        = 3;
    mIsothermMap["bet"]       = 3;
    mIsothermMap["quadratic"] = 3;
    mIsothermMap["dsl"]       = 4;
    mIsothermMap["dslf"]      = 6;
    }

IsothermModeler::~IsothermModeler()
    {
    }

std::shared_ptr<Isotherm>
IsothermModeler::fit(const std::string& isoname,
                     const std::vector<double>& x,
                     const std::vector<double>& y,
                     const std::vector<double>& guess)
    {
    if (mIsothermMap.count(isoname) == 0)
        throw IsothermException {__FILE__, __LINE__, "Invalid isotherm name."};

    if (x.size() != y.size())
        throw IsothermException {__FILE__, __LINE__, "Different data size."};

    if (x.empty())
        throw IsothermException {__FILE__, __LINE__, "Data size = 0."};

    IsothermFactory factory;

    int numParams = mIsothermMap[isoname];

    auto obj = [&isoname, &factory, &numParams, &x, &y](const Minimizer::PointType& p)
        {
        for (auto pi : p)
            if (pi < 1.0e-3 or pi > 1.0e3)
                return 1.0e30;

        std::vector<Any> params (numParams);
        for (int i = 0; i < numParams; ++i)
            params[i] = p[i];

        auto iso = factory.create(isoname, params);

        double error = 0.0;
        int size = x.size();
        for (int i = 0; i < size; ++i)
            {
            double dy = (iso->loading(x[i]) - y[i]) / (1.0 + y[i]);
            error += dy * dy;
            }

        return error / static_cast<double>(size);
        };

    Simplex simplex;

    std::vector<double> initPoint (numParams, 1.0);
    if (not guess.empty())
        initPoint = guess;

    simplex.setInitialPoint(initPoint).setFunction(obj);

    try {
        auto p = simplex.minimize().getMinimumPoint();

        std::vector<Any> params;
        for (auto x : p)
            params.push_back(x);

        auto iso = factory.create(isoname, params);

        double   yAcc = 0.0;
        double ysqAcc = 0.0;
        for (auto yi : y)
            {
              yAcc += yi;
            ysqAcc += yi * yi;
            }

        double stot = ysqAcc - yAcc / static_cast<double>(y.size());

        int xsize = x.size();
        double sres = 0.0;
        for (int i = 0; i < xsize; ++i)
            {
            double e = (iso->loading(x[i]) - y[i]);
            sres += e * e;
            }

        mRSquare = 1.0 - sres / stot;
        mError = simplex.getMinimumValue();

        return iso;
        }
    catch (MinimizerException& e)
        {
        std::string msg {"Isotherm fitting fails\n    Reason: "};
        msg += e.what();

        throw IsothermException {__FILE__, __LINE__, msg};
        }
    }

std::shared_ptr<Isotherm>
IsothermModeler::autofit(const std::vector<double>& x,
                         const std::vector<double>& y)
    {
    std::vector< std::shared_ptr<Isotherm> > isotherms;
    std::vector<double> errors;
    std::vector<double> rsquares;

    for (auto& info : mIsothermMap)
        {
        std::string isoname = info.first;
        double dimfactor = static_cast<double>(info.second); // second = # params

        isotherms.push_back(this->fit(isoname, x, y));
        errors.push_back(this->getError() * dimfactor);
        rsquares.push_back(this->getRSquare());
        }

    int index = std::min_element(errors.begin(), errors.end()) - errors.begin();

    // Check the R^2 of minimum error model.
    if (rsquares[index] < 0.99)
        {
        // Interpolator is the best.
        mRSquare = 0.0;
        mError = 0.0;
        std::vector<Any> args {x, y};
        return IsothermFactory {}.create("interpolator", args);
        }
    else
        {
        mRSquare = rsquares[index];
        mError = errors[index];
        return isotherms[index];
        }
    }

double
IsothermModeler::getError() const
    {
    return mError;
    }

double
IsothermModeler::getRSquare() const
    {
    return mRSquare;
    }
