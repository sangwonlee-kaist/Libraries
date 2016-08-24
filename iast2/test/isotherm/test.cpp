#include <iostream>
#include <iomanip>
#include <memory>
#include "../../isotherm.hpp"
#include "../../langmuir_isotherm.hpp"
#include "../../lf_isotherm.hpp"
#include "../../dsl_isotherm.hpp"
#include "../../dslf_isotherm.hpp"

int
main(int argc, char* argv[])
    {
    using namespace std;
    shared_ptr<Isotherm> iso = make_shared<LangmuirIsotherm>(1.0, 1.0);

    double p  = 2.0;
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

    return 0;
    }
