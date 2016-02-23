#include "iast.hpp"

#include <iostream>
#include <iomanip>

#define DEBUG(x) std::cout << #x << " = " << x << endl;

int
main()
try {
    iast_core myiast;
   
    interpolation_isotherm isotherm_array[] 
        {{"c1.txt"}, {"c2.txt"},{"c3.txt"}, {"c4.txt"}, {"c5.txt"}};  
       
    for (auto& iso : isotherm_array)
        {
        myiast.add_isotherm(iso);
        }
        
    myiast.set_temperature(300.0);
    myiast.set_composition({5.0, 4.0, 3.0, 2.0, 1.0});
    
    for (double P = 0.001; P < 1000.1; P *= pow(10, 1.0 / 10.0))
        {
        myiast.set_pressure(P);
        auto result = myiast.calculate();
        
        std::cout << P;
        for (auto ni : result.get_loadings())
            {
            std::cout << std::setw(15) << ni;
            }
        std::cout << std::setw(15) << result.get_cycle() << std::endl;
        }        
        
    return 0;
    }
catch (std::exception& e)
    {
    std::cout << e.what() << std::endl;
    }
