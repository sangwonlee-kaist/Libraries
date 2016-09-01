#include <iostream>
#include <iomanip>
#include <exception>
#include <vector>

#include "isotherm_utility.hpp"

int
main(int argc, char* argv[])
try {
    if (argc < 2)
        {
        std::cerr << "Usage: fitiso datafile." << std::endl;
        return 1;
        }

    IsothermModeler modeler;

    std::vector<double> x;
    std::vector<double> y;
    ::readTwoColumns(argv[1], x, y);
    auto iso = modeler.autofit(x, y);


    std::cout << iso->getInfoString() << std::endl;

    std::cout << std::setw(15) << "x" <<
                 std::setw(15) << "y" <<
                 std::setw(15) << "y fit" <<
                 std::endl;

    for (int i = 0; i < x.size(); ++i)
        std::cout << std::setw(15) << x[i] <<
                     std::setw(15) << y[i] <<
                     std::setw(15) << iso->loading(x[i]) <<
                     std::endl;

    return 0;
    }
catch (std::exception& e)
    {
    std::cerr << e.what() << std::endl;
    return 1;
    }
