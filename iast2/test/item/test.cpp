#include <iostream>
#include <functional>
#include <vector>

#include "../../isotherm_factory.hpp"
#include "../../isotherm_utility.hpp"
#include "../../linear_interpolator.hpp"

using namespace std;

int
main(int, char* [])
    {
    vector<double> x;
    vector<double> y;

    ::readTwoColumns("Q.dat", x, y);
    LinearInterpolator q;
    q.setData(x, y);

    function<double(double)> isoheat = q;

    IsothermFactory factory;

    ::readTwoColumns("n.dat", x, y);
    auto n = factory.create("interpolator", {x, y});

    auto iso = factory.create("item", {n, isoheat, 293.0, 303.0});

    for (double p = 0.0; p < 1000.0; p += 10.0)
        cout << p << "  " << iso->loading(p) << endl;

    iso = factory.create("item", {  n, isoheat, 293.0, 298.0});
    iso = factory.create("item", {iso, isoheat, 298.0, 303.0});

    for (double p = 0.0; p < 1000.0; p += 10.0)
        cout << p << "  " << iso->loading(p) << endl;

    return 0;
    }
