#pragma once

#include <exception>
#include <string>

class MinimizerException : public std::exception
    {
public:
    MinimizerException(const char* file, const int& line, const std::string& msg);

    virtual const char* what() const noexcept override;
private:
    std::string mMessage;
    };

MinimizerException::MinimizerException(
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
MinimizerException::what() const noexcept
    {
    return mMessage.c_str();
    }
