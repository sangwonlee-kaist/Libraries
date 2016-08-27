#include <iostream>
#include <iomanip>
#include <memory>
#include "../../isotherm.hpp"
#include "../../isotherm_factory.hpp"
#include "../../isotherm_exception.hpp"
#include "../../interpolator_isotherm.hpp"

int
main(int argc, char* argv[])
    {
    using namespace std;
    IsothermFactory factory;
    auto iso = factory.create("langmuir", {1.0, 1.0});

    double p  = 5.5;
    double dp = 1.e-6;

    cout << "1." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    iso = factory.create("lf", {1.0, 1.0, 2.0});

    cout << "2." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    iso = factory.create("dsl", {1.0, 1.0, 2.0, 2.0});

    cout << "3." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    iso = factory.create("dslf", {1.0, 1.0, 1.0, 2.0, 2.0, 2.0});

    cout << "4." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    vector<double> x {1.0, 2.0, 3.0, 4.0, 5.0};
    vector<double> y {0.5, 0.6666, 0.75, 0.8, 0.8333333};

    iso = factory.create("interpolator", {x, y});

    dynamic_pointer_cast<InterpolatorIsotherm>(iso)->pushBack(6.0, 0.857143);

    cout << "5." << endl;
    cout << iso->getInfoString() << endl;
    cout << iso->loading(p) << " = ";
    cout << (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;

    try {
        factory.create("what?", {0});
        }
    catch(IsothermException& e)
        {
        cout << e.what() << endl;
        }

    try {
        factory.create("lf", {1, 1, 3});
        }
    catch(IsothermException& e)
        {
        cout << e.what() << endl;
        }

    return 0;
    }
