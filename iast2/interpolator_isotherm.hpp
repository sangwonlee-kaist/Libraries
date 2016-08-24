#pragma once

#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cstddef>

#include "isotherm.hpp"
#include "linear_interpolator.hpp"

class InterpolatorIsotherm : public Isotherm
    {
public:
    InterpolatorIsotherm(const std::vector<double>& x,
                         const std::vector<double>& y);

    virtual ~InterpolatorIsotherm() = default;

    virtual double loading(double P) const override;
    virtual double spressure(double P) const override;

    virtual std::string getInfoString() const override;
private:
    mutable LinearInterpolator mLoading;
    std::vector<double> mSpressure;
    };

InterpolatorIsotherm::InterpolatorIsotherm(const std::vector<double>& x,
                                           const std::vector<double>& y)
    :
    mLoading {},
    mSpressure {}
    {
    if (x.size() != y.size())
        ;

    mLoading.setData(x, y);

    mSpressure.resize(x.size());
    mSpressure.front() = y.front();

    for (int i = 1; i < x.size(); ++i)
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
    const auto& x = mLoading.getXData();
    const auto& y = mLoading.getYData();

    if (P < x.front())
        return y.front() / x.front() * P;

    return mLoading(P);
    }

#include <iostream>
using namespace std;

double
InterpolatorIsotherm::spressure(double P) const
    {
    const auto& xdata = mLoading.getXData();
    const auto& ydata = mLoading.getYData();

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
    ss << "[Parameters] None";

    return ss.str();
    }
