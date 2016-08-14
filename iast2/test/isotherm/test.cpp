#include <iostream>
#include <memory>
#include "../../isotherm.hpp"

class MyIsotherm : public Isotherm
    {
public:
    virtual double loading(double P) override
        {
        return P;
        }

    virtual double spressure(double P) override
        {
        return P * P;
        }
    };

int
main(int argc, char* argv[])
    {
    using namespace std;
    shared_ptr<Isotherm> myIsotherm = make_shared<MyIsotherm>();

    cout << myIsotherm->loading(2.0) << endl;
    cout << myIsotherm->spressure(2.0) << endl;

    return 0;
    }
