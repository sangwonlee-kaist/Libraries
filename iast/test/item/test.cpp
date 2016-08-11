#include "../../../iast.hpp"
#include <iostream>
#include <fstream>
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
    InverseIsotherm in2 {n2, 298.0};
    for (double nn = 5.0; nn < 45.0; nn += 5.0)
        cout << in2(nn) << endl;

    cout << "Test 2-2" << endl;
    for (double nn = 0.1; nn < 10.0; nn += 0.1)
        cout << InverseIsotherm {n, 298.0}(nn) << "   " << nn << endl;

    cout << "Test 2-3" << endl;
    for (double nn = 0.1; nn < 10.0; nn += 0.1)
        cout << InverseIsotherm {n.get_loading(), 298.0}(nn) << "   " << nn << endl;

    cout << "Test 2-4" << endl;
    InverseIsotherm in {n.get_loading(), 298.0};
    for (double nn = 0.1; nn < 10.0; nn += 0.1)
        cout << in(nn) << endl;

    // Test 3 ====================================================
    ItemIsotherm itemed = ItemIsotherm {n, 293.0, "Q.dat"};

    cout << "Test 3-1" << endl;
    for (double p = 0.01; p <= 9.0; p += 0.01)
        cout << p << "    " << itemed.loading(303.0, p) << endl;

    double value = (itemed.spreading_pressure(303.0, 5.0 + 1.e-6) - itemed.spreading_pressure(303.0, 5.0)) / 1.e-6 * 5.0;
    cout << itemed.loading(303.0, 5.0) << "     " << value << endl;

    itemed.printPressureList();

    cout << "ZIF-8 CO2 test" << endl;

    DSLF_isotherm co2 {27.25, 2.122, 0.01533, 0.006895, 1.6, 0.3244};
    ItemIsotherm co2Item (co2, 298, "Qco2.dat");

    for (double t = 228.0; t < 298.0 + 0.1; t += 10.0)
        {
        ofstream ofs (to_string(t) + ".dat");
        for (double p = 0.1; p <= 25.0 + 0.1; p += 0.1)
            {
            ofs << p << "  " << co2Item.loading(t, p) << endl;
            }
        }

    return 0;
    }
