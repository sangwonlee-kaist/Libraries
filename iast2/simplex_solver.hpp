#pragma once

#include "solver.hpp"
#include "simplex.hpp"
#include <cmath>

class SimplexSolver : public Solver
    {
public:
    enum Option : int {NUM_REPEATS, TOL_X, TOL_F};

             SimplexSolver();
    virtual ~SimplexSolver() = default;

    virtual Solver& setFunctions(std::vector<FunctionType> functions) override;
    virtual Solver& setInitialPoint(PointType point) override;
    virtual Solver& setOption(int option, ValueType value) override;

    virtual Solver& solve() override;

    virtual PointType getRootPoint() const override;
    virtual int       getNumFunctionCalls() const override;
private:
    std::vector<FunctionType> mFunctions;
    PointType mInitialPoint;
    PointType mRootPoint;
    int mNumRepeats;
    ValueType mTolX;
    ValueType mTolF;
    int mNumFunctionCalls;
    };

SimplexSolver::SimplexSolver() :
    mFunctions {},
    mInitialPoint {},
    mRootPoint {},
    mNumRepeats {1},
    mTolX {1.0e-10},
    mTolF {1.0e-10},
    mNumFunctionCalls {}
    {

    }

Solver&
SimplexSolver::setFunctions(std::vector<FunctionType> functions)
    {
    mFunctions = functions;
    return *this;
    }

Solver&
SimplexSolver::setInitialPoint(PointType point)
    {
    mInitialPoint = point;
    return *this;
    }

Solver&
SimplexSolver::setOption(int option, ValueType value)
    {
    switch (option)
        {
        case Option::NUM_REPEATS:
            mNumRepeats = static_cast<int>(value + 1.0e-5);
            break;
        case Option::TOL_X:
            mTolX = value;
            break;
        case Option::TOL_F:
            mTolF = value;
            break;
        }

    return *this;
    }

Solver&
SimplexSolver::solve()
    {
    FunctionType objective = [this](const PointType& p)
        {
        double value = 0.0;
        for (auto& func : mFunctions)
            value += std::pow(func(p), 2);
        return value;
        };

    Simplex simplex;

    simplex.setFunction(objective).
        setInitialPoint(mInitialPoint).
        setOption(Simplex::Option::NUM_REPEATS, mNumRepeats).
        setOption(Simplex::Option::TOL_X, mTolX).
        setOption(Simplex::Option::TOL_F, mTolF);

    simplex.minimize();

    mRootPoint = simplex.getMinimumPoint();
    mNumFunctionCalls = simplex.getNumFunctionCalls();

    return *this;
    }

Solver::PointType
SimplexSolver::getRootPoint() const
    {
    return mRootPoint;
    }

int
SimplexSolver::getNumFunctionCalls() const
    {
    return mNumFunctionCalls;
    }
