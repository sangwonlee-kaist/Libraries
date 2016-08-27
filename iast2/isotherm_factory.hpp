#pragma once

#include <vector>
#include <memory>
#include <string>

#include "isotherm.hpp"
#include "any.hpp"

class IsothermFactory
    {
public:
     IsothermFactory() = default;
    ~IsothermFactory() = default;

    std::shared_ptr<Isotherm> create(const std::string& name, std::vector<Any> args) const;
    };
