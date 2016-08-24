#pragma once

#include <memory>
#include <string>

#include "simplex_solver.hpp"
#include "bisection_solver.hpp"
//#include "arma_solver.hpp"

class SolverFactory
    {
public:
             SolverFactory() = default;
    virtual ~SolverFactory() = default;

    virtual std::shared_ptr<Solver> create(std::string name);
    };

std::shared_ptr<Solver>
SolverFactory::create(std::string name)
    {
    if (name == "simplex")
        return std::make_shared<SimplexSolver>();
//    else if (name == "arma")
//        return std::make_shared<ArmaSolver>();
    else if (name == "bisection")
        return std::make_shared<BisectionSolver>();
    else
        throw SolverException {__FILE__, __LINE__, "Unsupported solver."};
    }
