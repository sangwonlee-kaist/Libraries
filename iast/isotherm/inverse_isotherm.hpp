#pragma once

class InverseIsotherm
    {
public:
    InverseIsotherm(const isotherm_base& isotherm,
                    const        double& temperature);

    //double pressure(double loading) const;
    double operator () (double loading) const;
private:
    const isotherm_base& mIsotherm;
    double               mTemperature;
    };

InverseIsotherm::InverseIsotherm(const isotherm_base& isotherm,
                                 const        double& temperature) :
    mIsotherm    {isotherm}, 
    mTemperature {temperature}
    {
    
    }

double
InverseIsotherm::operator () (double loading) const
    {
    // Simple Bisection Algorithm
    double xLow    = 0.0;
    double xHigh   = 0.0;
    double xCenter = 0.0;

    double fLow    = 0.0;
    double fHigh   = 0.0;
    double fCenter = 0.0;

    int iter = 0;
    int maxIter = 100;

    auto& iso = mIsotherm;
    auto& T   = mTemperature;

    xLow = 0.0;
    fLow = -loading;

    xHigh = 1.0;
    fHigh = iso.loading(T, xHigh) - loading;

    while (fHigh < 0.0)
        {
        xHigh *= 1.5;
        fHigh  = iso.loading(T, xHigh) - loading;

        if (xHigh > 10000.0)
            {
            throw std::runtime_error
                {"InverseIsotherm::inverseIsotherm(): Given uptake beyonds saturation loading."};
            }
        }

    maxIter = static_cast<int>(std::log2((xHigh - xLow)/ 1.e-6)) + 1;
    //std::cout << "max iter = " << maxIter << std::endl;
    for (iter = 0; iter < maxIter; ++iter)
        {
        xCenter = 0.5 * (xLow + xHigh);
        fCenter = iso.loading(T, xCenter) - loading;

        if (fCenter > 0.0)
            {
            xHigh = xCenter;
            fHigh = fCenter;
            }
        else if (fCenter < 0.0)
            {
            xLow = xCenter;
            fLow = fCenter;
            }
        }

    return xCenter;
    }
