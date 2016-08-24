#pragma once

#include <exception>
#include <string>

class IastException : public std::exception
    {
public:
    IastException(const char* file, const int& line, const std::string& msg);

    virtual const char* what() const noexcept override;
private:
    std::string mMessage;
    };

IastException::IastException(
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
IastException::what() const noexcept
    {
    return mMessage.c_str();
    }
