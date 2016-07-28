#pragma once

class interpolation_isotherm : public isotherm_base
    {
public:
    interpolation_isotherm() : p {}, q {}, pi {} {}
    interpolation_isotherm(const std::string& filename) :
        p {}, q {}, pi {}
        {
        std::ifstream ifs {filename};

        if (!ifs.good())
            {
            throw std::invalid_argument
                {"interpolation_isotherm() file does not exist."};
            }

        double P, Q;
        std::string single_line;
        while(std::getline(ifs, single_line))
            {
            std::stringstream ss;
            ss << single_line;
            ss >> P >> Q;
            p.push_back(P);
            q.push_back(Q);
            }

        ifs.close();

        size_t size {p.size()};
        pi.resize(size);
        // initial assumption near zero.
        pi[0] = q[0];
        // cache integration values.
        for (size_t i {1}; i < size; ++i)
            {
            double slope     {(q[i] - q[i - 1]) / (p[i] - p[i - 1])};
            double intercept {-slope * p[i] + q[i]};

            pi[i] = pi[i - 1] + slope * (p[i] - p[i - 1]) + intercept * std::log(p[i] / p[i - 1]);
            }
        }

    real_t
    loading(real_t T, real_t P) const
        override
        {
        if (P <= 0.0)
            {
            return 0.0;
            }

        size_t j = std::lower_bound(p.begin(), p.end(), P) - p.begin();

        if (j == 0)
            {
            return q[j] / p[j] * P;
            }

        if (j == q.size())
            {
            // Constant loading assumption.
            return q[j - 1];
            }

        double slope {(q[j] - q[j - 1]) / (p[j] - p[j - 1])};

        return slope * (P - p[j]) + q[j];
        }

    real_t
    spreading_pressure(real_t T, real_t P) const
        override
        {
        if (P <= 0.0)
            {
            return 0.0;
            }

        size_t j = std::lower_bound(p.begin(), p.end(), P) - p.begin();

        if (j == 0)
            {
            return q[j] / p[j] * P;
            }

        if (j == q.size())
            {
            return pi[j - 1] + q[j - 1] * std::log(P / p[j - 1]);
            }

        double slope     {(q[j] - q[j - 1]) / (p[j] - p[j - 1])};
        double intercept {-slope * p[j] + q[j]};

        return pi[j - 1] + slope * (P - p[j - 1]) + intercept * std::log(P / p[j - 1]);
        }

    real_t getMaxLoading() {return q.back();}
    real_t getMaxPressure() {return p.back();}
    real_t push_back(real_t _p, real_t _q)
        {
        p.push_back(_p);
        q.push_back(_q);

        pi.resize(p.size());

        if (p.size() == 1)
            {
            pi[0] = q[0];
            }
        else
            {
            int i = p.size() - 1;

            double slope     {(q[i] - q[i - 1]) / (p[i] - p[i - 1])};
            double intercept {-slope * p[i] + q[i]};

            pi[i] = pi[i - 1] + slope * (p[i] - p[i - 1]) + intercept * std::log(p[i] / p[i - 1]);
            }
        }
    bool empty() const
        {
        return p.empty();
        }
private:
    std::vector<double> p;
    std::vector<double> q;
    std::vector<double> pi;
    };




