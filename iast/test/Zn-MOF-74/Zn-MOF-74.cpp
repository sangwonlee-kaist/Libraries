#include "iast.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int
main(int argc, char* argv[])
try {
    DSLF_isotherm h2o {0.193274, 37.7745, 2.27e13, 3.19e27, 0.526687, 0.102092};
    DSLF_isotherm co2 {0.28502, 10.6243, 3.00e2, 1.83, 2.20e-1, 8.80e-1};

    iast_core binary;

    binary.add_isotherm(h2o);
    binary.add_isotherm(co2);

    binary.set_temperature(313.0);
    binary.set_pressure(0.15);

    for (double y {1.e-4}; y < 0.1; y *= pow(10.0, 0.1)) 
        {
        binary.set_composition({y, 1.0 - y});

        auto result = binary.calculate();

        cout << setw(15) << y;
        for (auto ld : result.get_loadings())
            {
            cout << setw(15) << ld;
            }
        cout << setw(15) << result.get_cycle() << endl;
        }
    
    return 0;
    }
catch (exception& e)
    {
    cout << e.what() << endl;
    }
