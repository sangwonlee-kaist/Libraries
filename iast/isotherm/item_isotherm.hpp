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
    double               mRefTemperature;
    const isotherm_base& mRefIsotherm;
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
    InverseIsotherm             mInverseIsotherm;

    void checkAndExpand(double T, double P) const;
public:
//    double inverseIsotherm(double loading) const;
    void   printPressureList() const;
    };

ItemIsotherm::ItemIsotherm(const isotherm_base& isothermBase,
                           const double&        refTemperature,
                           const std::string&   filename) :
    mRefTemperature  {refTemperature},
    mRefIsotherm     {isothermBase},
    mTargetIsotherms {},
    mIsostericHeat   {filename},
    mInverseIsotherm {mRefIsotherm, mRefTemperature}
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
    auto& invIso = mInverseIsotherm;
    auto& PList  = mPressureList;
    auto& Q      = mIsostericHeat;

    double tiny  = 1.e-10;
    double R     = 0.008314469;
    double dbeta = (1.0 / T - 1.0 / Tref) / R;

    if (iso.empty())
        {
        // push first data.
        if (PList.empty())
            PList.push_back(invIso(0.1));

        double Ptar = PList[0] * std::exp(-dbeta * Q(0.1));
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
                PList.push_back(invIso(nn));
            }

        double Ptar = PList[index] * std::exp(-dbeta * Q(n));
        iso.push_back(Ptar, n);
        }
    }

void
ItemIsotherm::printPressureList() const
    {
    for (auto& p : mPressureList)
        std::cout << p << std::endl;
    }
