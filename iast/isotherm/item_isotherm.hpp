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
    RealType            loading(RealType T, RealType P) override;
    RealType spreading_pressure(RealType T, RealType P) override;
private:
    double               mRefTemperature;
    const isotherm_base& mRefIsotherm;
    IsothermContainer    mTargetIsotherms;
    Interpolator         mIsostericHeat;
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
ItemIsotherm::loading(RealType T, RealType P)
    {
    std::string Tstr = std::to_string(T);

    if (Tstr == std::to_string(mRefTemperature))
     //   return mRefrLoading(T, P);

    if (mTargetIsotherms.count(Tstr) == 0)
        {
        mTargetIsotherms[Tstr] = interpolation_isotherm {};
        }

    interpolation_isotherm& iso = mTargetIsotherms[Tstr];

    if (P > iso.getMaxPressure())
        {
        // Expand data.

        }

    return 0.0;
    }

ItemIsotherm::RealType
ItemIsotherm::spreading_pressure(RealType T, RealType P)
    {
    std::string Tstr = std::to_string(T);

    if (Tstr == std::to_string(mRefTemperature))
    //    return mRefSpreadingPressure(T, P);

    return 0.0;
    }
