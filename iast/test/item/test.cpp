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

    //cout << test.getMaxPressure() << endl;

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
    cout << "Test 2-1" << endl;
    interpolation_isotherm n2 {"h2o.dat"};
    for (double nn = 5.0; nn < 45.0; nn += 5.0)
        cout << ItemIsotherm {n2, 293.0, "Q.dat"}.inverseIsotherm(nn) << ", " << nn << endl;

    cout << "Test 2-2" << endl;
    for (double nn = 0.1; nn < 10.0; nn += 0.1)
        cout << ItemIsotherm {n, 293.0, "Q.dat"}.inverseIsotherm(nn) << "   " << nn << endl;

    // Test 3 ====================================================
    ItemIsotherm itemed = ItemIsotherm {n, 293.0, "Q.dat"};

    cout << "Test 3-1" << endl;
    for (double p = 0.01; p <= 9.0; p += 0.01)
        cout << p << "    " << itemed.loading(303.0, p) << endl;

    double value = (itemed.spreading_pressure(303.0, 5.0 + 1.e-6) - itemed.spreading_pressure(303.0, 5.0)) / 1.e-6 * 5.0;
    cout << itemed.loading(303.0, 5.0) << "     " << value << endl;

    itemed.printPressureList();

    return 0;
    }
