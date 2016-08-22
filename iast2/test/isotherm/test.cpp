#include <iostream>
#include <iomanip>
#include <memory>
#include "../../isotherm.hpp"
#include "../../langmuir_isotherm.hpp"

int
main(int argc, char* argv[])
    {
    using namespace std;
    shared_ptr<Isotherm> iso = make_shared<LangmuirIsotherm>(1.0, 1.0);

    double p  = 2.0;
    double dp = 1.e-6;

    cout << setw(15) << "1. Langmuir isotherm" << endl;
    cout << setw(15) << iso->loading(p) << endl;
    cout << setw(15) <<
        (iso->spressure(p + dp) - iso->spressure(p)) / dp * p << endl;


    return 0;
    }
