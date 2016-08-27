#include "interpolator_isotherm.hpp"

#include <cmath>
#include <sstream>
#include <algorithm>
#include <cstddef>

#include "linear_interpolator.hpp"

InterpolatorIsotherm::InterpolatorIsotherm(const std::vector<double>& x,
                                           const std::vector<double>& y)
    :
    mLoading {std::make_shared<LinearInterpolator>()},
    mSpressure {}
    {
    //if (x.size() != y.size())
    //    ;

    mLoading->setData(x, y);

    mSpressure.resize(x.size());
    mSpressure.front() = y.front();

    int maxi = static_cast<int>(x.size());
    for (int i = 1; i < maxi; ++i)
        {
        // y = a * x + b;
        double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        double b = -a * x[i] + y[i];

        // integrate(y / x) = a * x + log(x) + C
        mSpressure[i] = mSpressure[i - 1] +
            a * (x[i] - x[i - 1]) + b * std::log(x[i] / x[i - 1]);
        }
    }

double
InterpolatorIsotherm::loading(double P) const
    {
    const auto& x = mLoading->getXData();
    const auto& y = mLoading->getYData();

    if (P < x.front())
        return y.front() / x.front() * P;

    return (*mLoading)(P);
    }

double
InterpolatorIsotherm::spressure(double P) const
    {
    const auto& xdata = mLoading->getXData();
    const auto& ydata = mLoading->getYData();

    std::size_t i =
        std::lower_bound(xdata.begin(), xdata.end(), P) - xdata.begin();

    if (i == 0 and P < xdata.front())
        {
        const double x = xdata.front();
        const double y = ydata.front();

        return y / x * P;
        }

    if (i == mSpressure.size())
        {
        const double x = xdata.back();
        const double y = ydata.back();

        return mSpressure.back() + y * std::log(P / x);
        }

    const auto& x = xdata;
    const auto& y = ydata;

    double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    double b = -a * x[i] + y[i];

    return mSpressure[i - 1] + a * (P - x[i - 1]) + b * std::log(P / x[i - 1]);
    }

std::string
InterpolatorIsotherm::getInfoString() const
    {
    std::stringstream ss;

    ss << "[Interpolator Isotherm]\n";
    ss << "[Parameters] front (x, y) = (" <<
          mLoading->getXData().front() << ", " <<
          mLoading->getYData().front() << ") " <<
          ", back (x, y) = (" <<
          mLoading->getXData().back() << ", " <<
          mLoading->getYData().back() << ") ";

    return ss.str();
    }

void
InterpolatorIsotherm::pushBack(double p, double n)
    {
    mLoading->pushBack(p, n);

    const auto& x = mLoading->getXData();
    const auto& y = mLoading->getYData();

    // last index.
    int i = x.size() - 1;

    // y = a * x + b;
    double a = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    double b = -a * x[i] + y[i];

    // integrate(y / x) = a * x + log(x) + C
    mSpressure.push_back(mSpressure.back() +
        a * (x[i] - x[i - 1]) + b * std::log(x[i] / x[i - 1]));
    }
