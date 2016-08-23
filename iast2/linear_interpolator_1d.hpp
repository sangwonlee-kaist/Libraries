#pragma once

#include "interpolator.hpp"

class LinearInterpolator1D : Interpolator1D
    {
public:
             LinearInterpolator1D();
    virtual ~LinearInterpolator1D() = default;

    virtual Interpolator1D& setData(std::vector<DataPointType> data) override;
    virtual double operator () (const PointType& x) override;
private:
    std::vector<DataPointType> mData;
    };

LinearInterpolator1D::LinearInterpolator1D() :
    x {}, y {}
    {

    }

Interpolator1D&
LinearInterpolator1D::setData(std::vector<DataPointType> data)
    {
    mData = data;
    }


