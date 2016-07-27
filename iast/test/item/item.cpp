#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>

class Interpolation
    {
public:
    Interpolation(const std::string& filename);
    ~Interpolation();

    double operator () (double x);
private:
    std::vector<double> xs;
    std::vector<double> ys;
    };

Interpolation::Interpolation(const std::string& filename)
    {
    std::ifstream ifs {filename};

    double x, y;
    std::string buffer;
    while(std::getline(ifs, buffer))
        {
        std::stringstream ss;
        ss << buffer;
        ss >> x >> y;
        xs.push_back(x);
        ys.push_back(y);
        }

    ifs.close();
    }

Interpolation::~Interpolation()
    {

    }

double
Interpolation::operator () (double x)
    {
    std::size_t i = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();

    if (i == 0)
        return ys[i] / xs[i] * x;

    if (i == ys.size())
        return ys[i - 1];

    double slope = (ys[i] - ys[i - 1]) /
                   (xs[i] - xs[i - 1]);

    return slope * (x - xs[i]) + ys[i];
    }

int
main(int, char* [])
    {
    Interpolation q ("Q.dat");
    Interpolation n ("n.dat");

    const double R = 0.008314469;

    for (double dT = 10.0; dT < 180.1; dT += 10.0)
        {
        double T = 293.0 + dT;
        std::ofstream ofs (std::to_string(T) + ".dat");

        // delta beta.
        double db = 1.0 / (293.0 + dT) / R - 1.0 / 293.0 / R;

        for (double P = 0.0; P < 1000.0; P += 10.0)
            {
            double oldN, newN;

            oldN = n(P);
            for (int i = 0; i < 100; ++i)
                {
                newN = n(P * std::exp(q(oldN) * db));

                if (std::abs(1.0 - oldN / newN) < 1.e-5)
                    break;

                oldN = newN;
                }

            ofs << std::setw(15) << P << std::setw(15) << newN << std::endl;
            }
        } // dT

    return 0;
    }
