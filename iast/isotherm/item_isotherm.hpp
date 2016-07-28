#pragma once

class Interpolator
    {
public:
    Interpolator();
    Interpolator(const std::string& filename);
    Interpolator(std::ifstream& ifs);
    ~Interpolator();

    void setData(const std::string& filename);
    void setData(std::ifstream& ifs);
    double operator () (double x) const;
private:
    std::vector<double> xs;
    std::vector<double> ys;
    };

Interpolator::Interpolator()
    {

    }

Interpolator::Interpolator(const std::string& filename)
    {
    this->setData(filename);
    }

Interpolator::Interpolator(std::ifstream& ifs)
    {
    this->setData(ifs);
    }

Interpolator::~Interpolator()
    {

    }

void
Interpolator::setData(const std::string& filename)
    {
    std::ifstream ifs {filename};

    if (not ifs)
        throw std::invalid_argument {"Interpolator: invalid data file."};

    this->setData(ifs);
    }

void
Interpolator::setData(std::ifstream& ifs)
    {
    xs.resize(0);
    ys.resize(0);

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
    }

double
Interpolator::operator () (double x) const
    {
    std::size_t i = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();

    if (i == 0)
        return ys[i];

    if (i == ys.size())
        return ys[i - 1];

    double slope = (ys[i] - ys[i - 1]) /
                   (xs[i] - xs[i - 1]);

    return slope * (x - xs[i]) + ys[i];
    }



// Add ITEM functionality.
// A Decolator class.
class ItemIsotherm : public isotherm_base
    {
public:
    using RealType          = isotherm_base::real_t;
    using IsothermContainer = std::map<std::string, interpolation_isotherm>;

    ItemIsotherm(const isotherm_base& isothermBase,
                 const double&        refTemperature,
                 const std::string&   filename);
    RealType            loading(RealType T, RealType P) const override;
    RealType spreading_pressure(RealType T, RealType P) const override;
private:
    double                    mRefTemperature;
    const isotherm_base&      mRefIsotherm;
    // PressureList
    // P | Uptake
    // --+-------
    // P0|  0.1
    // P1|  0.2
    // P3|  0.3
    // ...
    // Here, Pi is PressureList!
    // Which is obtained using reference isotherm.
    // We can access to element using
    mutable std::vector<double> mPressureList;
    mutable IsothermContainer   mTargetIsotherms;
    Interpolator                mIsostericHeat;

    void checkAndExpand(double T, double P) const;
public:
    double inverseIsotherm(double loading) const;
    void   printPressureList() const;
    };

ItemIsotherm::ItemIsotherm(const isotherm_base& isothermBase,
                           const double&        refTemperature,
                           const std::string&   filename) :
    mRefTemperature  {refTemperature},
    mRefIsotherm     {isothermBase},
    mTargetIsotherms {},
    mIsostericHeat   {filename}
    {

    }

ItemIsotherm::RealType
ItemIsotherm::loading(RealType T, RealType P) const
    {
    auto& Tref = mRefTemperature;
    auto& iso  = mRefIsotherm;
    std::string Tstr = std::to_string(T);

    if (Tstr == std::to_string(Tref))
       return iso.loading(T, P);

    this->checkAndExpand(T, P);

    return mTargetIsotherms[Tstr].loading(T, P);
    }

ItemIsotherm::RealType
ItemIsotherm::spreading_pressure(RealType T, RealType P) const
    {
    auto& Tref = mRefTemperature;
    auto& iso  = mRefIsotherm;
    std::string Tstr = std::to_string(T);

    if (Tstr == std::to_string(Tref))
        return iso.spreading_pressure(T, P);

    this->checkAndExpand(T, P);

    return mTargetIsotherms[Tstr].spreading_pressure(T, P);
    }

void
ItemIsotherm::checkAndExpand(double T, double P) const
    {
    auto& Tref = mRefTemperature;
    std::string Tstr = std::to_string(T);
    // Check temperature data exist.
    if (mTargetIsotherms.count(Tstr) == 0)
        mTargetIsotherms[Tstr] = interpolation_isotherm {};

    interpolation_isotherm& iso = mTargetIsotherms[Tstr];
    auto& PList = mPressureList;
    auto& Q     = mIsostericHeat;

    double tiny  = 1.e-10;
    double R     = 0.008314469;
    double dbeta = (1.0 / T - 1.0 / Tref) / R;

    if (iso.empty())
        {
        // push first data.
        if (PList.empty())
            PList.push_back(this->inverseIsotherm(0.1));

        double Ptar = PList[0] * std::exp(-dbeta * Q(0.1));
        //double Ptar = (this->inverseIsotherm(0.1)) * std::exp(-dbeta * Q(0.1));
        iso.push_back(Ptar, 0.1);
        }

    while (P > iso.getMaxPressure())
        {
        // Expand data.
        double n = iso.getMaxLoading() + 0.1;
        int index = static_cast<int>(n / 0.1 - tiny);

        if (PList.size() <= index)
            {
            for (double nn = (PList.size() + 1) * 0.1; nn <= n + tiny; nn += 0.1)
                PList.push_back(this->inverseIsotherm(nn));
            }

        double Ptar = PList[index] * std::exp(-dbeta * Q(n));
        //double Ptar = (this->inverseIsotherm(n)) * std::exp(-dbeta * Q(n));
        iso.push_back(Ptar, n);
        }
    }

double
ItemIsotherm::inverseIsotherm(double loading) const
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

    auto& iso = mRefIsotherm;
    auto& T   = mRefTemperature;

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
                {"ItemIsotherm::inverseIsotherm(): Given uptake beyonds saturation loading."};
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

// It would be nice if we implement Invertible decolator.
//double
//ItemIsotherm::inverseIsotherm(double loading) const
//    {
//    // Ridders' Method.
//    double P1 = 0.0;
//    double P2 = 0.0;
//    double P3 = 0.0;
//    double P4 = 0.0;
//    double N1 = 0.0;
//    double N2 = 0.0;
//    double N3 = 0.0;
//    double N4 = 0.0;
//    double bestP = 0.0;
//    int iter = 0;
//    int maxIter = 100;
//
//    auto& iso = mRefIsotherm;
//    auto& T   = mRefTemperature;
//
//    bestP = -1.0e30;
//
//    P1 = 0.0;
//    N1 = -loading;
//
//    // Find proper P2.
//    P2 = 1.0;
//    N2 = iso.loading(T, P2) - loading;
//
//     while (N2 < 0.0)
//        {
//        P2 *= 1.5;
//        N2  = iso.loading(T, P2) - loading;
//
//        if (P2 > 10000.0)
//            {
//            throw std::runtime_error
//                {"ItemIsotherm::inverseIsotherm(): Given uptake beyonds saturation loading."};
//            }
//        }
//
//    //std::cout << "P2 = " << P2 << std::endl;
//    //std::cout << "N2 = " << N2 << std::endl;
//
//    // P1 >> P2 >> P3 >> P4 ...
//    for (iter = 0; iter < maxIter; ++iter)
//        {
//        P3 = 0.5 * (P1 + P2);
//        N3 = iso.loading(T, P3) - loading;
//
//        double sign = N1 > N2 ? 1.0 : - 1.0;
//        double factor = N3 * N3 - N1 * N2;
//        //std::cout << "factor = " << factor << std::endl;
//        if (factor < 1.e-10)
//            break;
//
//        P4 = P3 + (P3 - P1) * sign * N3 / std::sqrt(factor);
//        //std::cout << iter << ", " << P4 << std::endl;
//
//        if (std::abs(1.0 - bestP / P4) < 1.e-10)
//            break;
//
//        N4 = iso.loading(T, P4) - loading;
//        if (std::abs(N4) < 1.e-7)
//            break;
//
//        // Update for next iteration.
//        if (std::copysign(N3, N4) != N3)
//            {
//            P1 = P3;
//            N1 = N3;
//            P2 = P4;
//            N2 = N4;
//            }
//        else if (std::copysign(N1, N4) != N1)
//            {
//            P2 = P4;
//            N2 = N4;
//            }
//        else if (std::copysign(N2, N4) != N2)
//            {
//            P1 = P4;
//            N1 = P4;
//            }
//        else
//            {
//            throw std::runtime_error
//                {"ItemIsotherm::inverseIsotherm(): Weird!"};
//            }
//
//        bestP = P4;
//
//        if (std::abs(1.0 - P1 / P2) < 1.e-10)
//            break;
//
//        //std::cout << P4 << std::endl;
//        }
//
//    if (iter == maxIter)
//        {
//        throw std::runtime_error
//            {"ItemIsotherm::inverseIsotherm: inversing fails."};
//        }
//
//    return P4;
//    }

void
ItemIsotherm::printPressureList() const
    {
    for (auto& p : mPressureList)
        std::cout << p << std::endl;
    }
