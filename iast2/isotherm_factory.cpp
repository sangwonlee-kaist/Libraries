#include "isotherm_factory.hpp"

#include "isotherm_exception.hpp"
#include "langmuir_isotherm.hpp"
#include "lf_isotherm.hpp"
#include "dsl_isotherm.hpp"
#include "dslf_isotherm.hpp"
#include "interpolator_isotherm.hpp"

std::shared_ptr<Isotherm>
IsothermFactory::create(const std::string& name, std::vector<Any> args) const
    {
    try {
        if (name == "langmuir")
            {
            double q1 = args[0].getAs<double>();
            double k1 = args[1].getAs<double>();

            return std::make_shared<LangmuirIsotherm>(q1, k1);
            }
        if (name == "lf")
            {
            double q1 = args[0].getAs<double>();
            double k1 = args[1].getAs<double>();
            double n1 = args[2].getAs<double>();

            return std::make_shared<LfIsotherm>(q1, k1, n1);
            }
        if (name == "dsl")
            {
            double q1 = args[0].getAs<double>();
            double k1 = args[1].getAs<double>();
            double q2 = args[2].getAs<double>();
            double k2 = args[3].getAs<double>();

            return std::make_shared<DslIsotherm>(q1, k1, q2, k2);
            }
        if (name == "dslf")
            {
            double q1 = args[0].getAs<double>();
            double k1 = args[1].getAs<double>();
            double n1 = args[2].getAs<double>();
            double q2 = args[3].getAs<double>();
            double k2 = args[4].getAs<double>();
            double n2 = args[5].getAs<double>();

            return std::make_shared<DslfIsotherm>(q1, k1, n1, q2, k2, n2);
            }
        if (name == "interpolator")
            {
            std::vector<double> x = args[0].getAs<std::vector<double>>();
            std::vector<double> y = args[1].getAs<std::vector<double>>();

            return std::make_shared<InterpolatorIsotherm>(x, y);
            }
        else
            {
            throw IsothermException {__FILE__, __LINE__, name + ": Unsupported isotherm."};
            }
        }
    catch (AnyException& e)
        {
        throw IsothermException {__FILE__, __LINE__, "Invalid arguments type for " + name};
        }
    }
