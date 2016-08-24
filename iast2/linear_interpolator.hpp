#pragma once

#include <algorithm>
#include <tuple>
#include <cstddef>

#include "interpolator.hpp"

class LinearInterpolator : public Interpolator
    {
public:
             LinearInterpolator() = default;
    virtual ~LinearInterpolator() = default;

    virtual std::vector<double>& getXData() override;
    virtual std::vector<double>& getYData() override;
    virtual double at(double x) const override;
    virtual double operator () (double x) const override;

    virtual Interpolator& setData(const std::vector<double>& x,
                                  const std::vector<double>& y) override;
    virtual Interpolator& pushBack(double x, double y) override;
    
    virtual Interpolator& setOption(int option, double value) override;
private:
    std::vector<double> xs;
    std::vector<double> ys;
    };

double
LinearInterpolator::at(double x) const
    {
    if (xs.empty())
        ;

    std::size_t i = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();

    if (i == 0)
        return ys.front();

    if (i == xs.size())
        return ys.back();

    double slope = (ys[i] - ys[i - 1]) /
                   (xs[i] - xs[i - 1]);

    return slope * (x - xs[i]) + ys[i];
    }

std::vector<double>&
LinearInterpolator::getXData()
    {
    return xs;
    }

std::vector<double>&
LinearInterpolator::getYData()
    {
    return ys;
    }

double
LinearInterpolator::operator () (double x) const
    {
    return this->at(x);
    }

Interpolator&
LinearInterpolator::setData(const std::vector<double>& x,
                            const std::vector<double>& y)
    {
    xs = x;
    ys = y;

    return *this;
    }

Interpolator&
LinearInterpolator::pushBack(double x, double y)
    {
    xs.push_back(x);
    ys.push_back(y);

    return *this;
    }

Interpolator&
LinearInterpolator::setOption(int option, double value)
    {
    // throw
    return *this;
    }
