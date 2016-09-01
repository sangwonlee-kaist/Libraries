#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <memory>
#include <map>

#include "isotherm.hpp"

void readTwoColumns(const std::string& filename,
                    std::vector<double>& x,
                    std::vector<double>& y);

double inverseIsotherm(Isotherm& isotherm, double n);

class IsothermModeler final
    {
public:
     IsothermModeler();
    ~IsothermModeler();

    std::shared_ptr<Isotherm> fit(const std::string& isoname,
                                  const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& guess =
                                        std::vector<double> {});

    std::shared_ptr<Isotherm> autofit(
                                  const std::vector<double>& x,
                                  const std::vector<double>& y);

    double getError() const;
    double getRSquare() const;
private:
    std::map<std::string, int> mIsothermMap;
    double mError;
    double mRSquare;
    };
