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
    double operator () (double x);
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
Interpolator::operator () (double x)
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
    mutable IsothermContainer mTargetIsotherms;
    Interpolator              mIsostericHeat;

    void checkAndExpand(const std::string& Tstr, double P) const;
public:
    double inverseIsotherm(double loading) const;
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
    std::string Tstr = std::to_string(T);

    if (Tstr == std::to_string(mRefTemperature))
       return mRefIsotherm.loading(T, P);

    this->checkAndExpand(Tstr, P);

    return 0.0;
    }

ItemIsotherm::RealType
ItemIsotherm::spreading_pressure(RealType T, RealType P) const
    {
    std::string Tstr = std::to_string(T);

    if (Tstr == std::to_string(mRefTemperature))
        return mRefIsotherm.spreading_pressure(T, P);

    this->checkAndExpand(Tstr, P);

    return 0.0;
    }

void
ItemIsotherm::checkAndExpand(const std::string& Tstr, double P) const
    {
    // Check temperature data exist.
    if (mTargetIsotherms.count(Tstr) == 0)
        {
        mTargetIsotherms[Tstr] = interpolation_isotherm {};
        }

    interpolation_isotherm& iso = mTargetIsotherms[Tstr];

    if (iso.empty())
        {
        // push first data.

        }

    while (P > iso.getMaxPressure())
        {
        // Expand data.
        }

    }

// It would be nice if we implement Invertible decolator.
double
ItemIsotherm::inverseIsotherm(double loading) const
    {
    // Secant Method.
    double oldP = 1.0;
    double newP = 0.0;
    double oldN = 0.0;
    double newN = 0.0;
    int iter = 0;
    int maxIter = 100;

    double loadingGuess = mRefIsotherm.loading(mRefTemperature, oldP);
    while (loadingGuess < loading)
        {
        oldP *= 1.5;
        loadingGuess = mRefIsotherm.loading(mRefTemperature, oldP);
        }

    oldN = 0.0;
    newN = loadingGuess;

    // oldP >> newP >> nextP ...
    for (iter = 0; iter < maxIter; ++iter)
        {
        std::cout << oldP << std::endl;
        
        double slope = (newN - oldN) / (newP - oldP);
        double nextP = newP - newN / slope;
        //newP = oldP / mRefIsotherm.loading(mRefTemperature, oldP) * loading;
        if (std::abs(1.0 - nextP / newP) < 1.e-4)
            break;

        oldP = newP;
        oldN = newN;
        newP = nextP;
        newN = mRefIsotherm.loading(mRefTemperature, nextP);
        }

    if (iter == maxIter)
        {
        loadingGuess = mRefIsotherm.loading(mRefTemperature, newP);
        if (std::abs(1.0 - loadingGuess / loading) > 0.1)
            throw std::runtime_error
                {"ItemIsotherm::inverseIsotherm: inversing fails."};
        }

    return newP;
    }
