#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <tuple>

#include "../../iast.hpp"
#include "../../langmuir_isotherm.hpp"
#include "../../lf_isotherm.hpp"
#include "../../interpolator_isotherm.hpp"

using namespace std;

void readTwoColumn(string name, vector<double>& x, vector<double>& y)
    {
    ifstream ifs {name};

    x.clear();
    y.clear();

    while (ifs)
        {
        double xx, yy;

        if (ifs >> xx >> yy)
            {
            x.push_back(xx);
            y.push_back(yy);
            }
        }
    }

int
main(int argc, char* argv[])
    {
    try {
        Iast::IsothermVector isotherms;
        isotherms.resize(4);

        isotherms[0] = make_shared<LangmuirIsotherm>(32.9301, 0.00924395); // n2
        isotherms[1] = make_shared<LfIsotherm>(10.4655, 1.90362, 1.1976);  // co2

        vector<double> x, y;
        readTwoColumn("h2o.dat", x, y);
        isotherms[2] = make_shared<InterpolatorIsotherm>(x, y);            // h2o
        isotherms[3] = make_shared<LangmuirIsotherm>(11.6083, 0.0220011);  // o2

        for (const auto& iso : isotherms)
            cout << iso->getInfoString() << endl;

        Iast iast;
        iast.setIsotherms(isotherms);

        vector<double> gasComposition {0.799, 0.15, 0.001, 0.05};

        double uptake;
        vector<double> composition;

        iast.calculate(Iast::Mode::FIX_PY, 1.0, gasComposition);
        tie(uptake, composition) = iast.getResult();

        cout << "Uptake = " << uptake << endl;
        cout << "(";
        for (auto& comp : composition)
            cout << comp << ", ";
        cout << ")";
        }
    catch (IastException& e)
        {
        cout << e.what() << endl;
        }

    return 0;
    }
