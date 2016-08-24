#include <iostream>
#include <iomanip>
#include <memory>
#include "../../isotherm.hpp"
#include "../../langmuir_isotherm.hpp"
#include "../../lf_isotherm.hpp"
#include "../../dsl_isotherm.hpp"
#include "../../dslf_isotherm.hpp"
#include "../../interpolator_isotherm.hpp"

int
main(int argc, char* argv[])
    {
    using namespace std;
    shared_ptr<Isotherm> iso = make_shared<LangmuirIsotherm>(1.0, 1.0);

    double p  = 10.6;
    double dp = 1.e-6;

    cout << "1." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    iso = make_shared<LfIsotherm>(1.0, 1.0, 2.0);

    cout << "2." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    iso = make_shared<DslIsotherm>(1.0, 1.0, 2.0, 2.0);

    cout << "3." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    iso = make_shared<DslfIsotherm>(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);

    cout << "4." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    std::vector<double> x {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y {0.5, 0.6666, 0.75, 0.8, 0.8333333};

    iso = make_shared<InterpolatorIsotherm>(x, y);

    cout << "5." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    return 0;
    }
