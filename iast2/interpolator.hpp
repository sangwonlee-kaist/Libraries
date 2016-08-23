#pragma once

#include <array>
#include <vector>
#include <cstddef>

template <std::size_t D>
class Interpolator
    {
public:
    using ValueType = double;
    using PointType = std::array<ValueType, D>;
    using DataPointType = std::array<ValueType, D + 1>;

             Interpolator() = default;
    virtual ~Interpolator() = default;

    virtual Interpolator<D>& setOption(int option, ValueType value) = 0;
    virtual Interpolator<D>& setData(std::vector<DataPointType> data) = 0;
    virtual double operator () (const PointType& x) = 0;
    };

using Interpolator1D = Interpolator<1>;
using Interpolator2D = Interpolator<2>;
using Interpolator3D = Interpolator<3>;
