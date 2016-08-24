#pragma once

#include <exception>
#include <string>

class SolverException : public std::exception
    {
public:
    SolverException(const char* file, const int& line, const std::string& msg);

    virtual const char* what() const noexcept override;
private:
    std::string mMessage;
    };

SolverException::SolverException(
    const char* file,
    const int&  line,
    const std::string& msg) : mMessage {}
    {
    mMessage += "File: ";
    mMessage += file;
    mMessage += ", Line: ";
    mMessage += std::to_string(line);
    mMessage += ": " + msg;
    }

const char*
SolverException::what() const noexcept
    {
    return mMessage.c_str();
    }


class NoRootException : public SolverException
    {
    using SolverException::SolverException;
    };
