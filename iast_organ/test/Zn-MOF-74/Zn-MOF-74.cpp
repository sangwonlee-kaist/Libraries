#include "iast.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int
main(int argc, char* argv[])
try {
    interpolation_isotherm h2o {"h2o.txt"};
    interpolation_isotherm co2 {"co2.txt"};

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
        cout << endl;
        }
    
    return 0;
    }
catch (exception& e)
    {
    cout << e.what() << endl;
    }
