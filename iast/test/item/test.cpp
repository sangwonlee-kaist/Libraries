#include "../../../iast.hpp"
#include <iostream>
using namespace std;

int
main()
    {
    // Test 1 ====================================================
    // interpolation_isotherm interface.
    interpolation_isotherm test;
    interpolation_isotherm n {"n.dat"};

    cout << n.getMaxPressure() << endl;
    cout << n.getMaxLoading() << endl;
    n.push_back(10010.0, 15.00);
    cout << n.getMaxPressure() << endl;
    cout << n.getMaxLoading() << endl;

    double  P = 10005.0;
    double dP = 1.e-5;
    cout << n.loading(0.0, P) << endl;
    cout << (n.spreading_pressure(0.0, P + dP) -
             n.spreading_pressure(0.0, P     )) / dP * P << endl;


    // Test 2 ====================================================
    ItemIsotherm {n, 293.0, "Q.dat"};

    return 0;
    }
