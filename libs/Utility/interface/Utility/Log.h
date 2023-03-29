#pragma once

#include <fmt/core.h>

#include <iomanip>
#include <iostream>

namespace plasmatic::Log {

enum class Level { Debug = 0, Info, Warn, Error };

namespace detail {
// Logger singleton class (do not use this directly!)
// See https://stackoverflow.com/a/1008289 for more info
class LoggerSingleton {
  public:
    LoggerSingleton(LoggerSingleton const &) = delete;
    void operator=(LoggerSingleton const &) = delete;

    static LoggerSingleton &getInstance() { return getInstanceImpl(); }

    void SetVerbosityLevel(const Level &level) { _level = level; }

    template <typename... Args> void Info(fmt::format_string<Args...> fmt, Args &&...args) {
        if (_level > Level::Info) {
            return;
        }

        std::time_t t = std::time(nullptr);
        std::cout << std::put_time(std::localtime(&t), "%F %T %Z") << " \u001b[32m[Info]\u001b[0m "
                  << fmt::format(fmt, std::forward<Args>(args)...) << std::endl;
    }

    template <typename... Args> void Warn(fmt::format_string<Args...> fmt, Args &&...args) {
        if (_level > Level::Warn) {
            return;
        }

        std::time_t t = std::time(nullptr);
        std::cout << std::put_time(std::localtime(&t), "%F %T %Z") << " \u001b[33m[Warn]\u001b[0m "
                  << fmt::format(fmt, std::forward<Args>(args)...) << std::endl;
    }

    template <typename... Args> void Error(fmt::format_string<Args...> fmt, Args &&...args) {
        if (_level > Level::Error) {
            return;
        }

        std::time_t t = std::time(nullptr);
        std::cerr << std::put_time(std::localtime(&t), "%F %T %Z") << " \u001b[31m[Error]\u001b[0m "
                  << fmt::format(fmt, std::forward<Args>(args)...) << std::endl;
    }

    template <typename... Args> void Debug(fmt::format_string<Args...> fmt, Args &&...args) {
        if (_level > Level::Debug) {
            return;
        }

        std::time_t t = std::time(nullptr);
        std::cout << std::put_time(std::localtime(&t), "%F %T %Z") << " \u001b[36m[Debug]\u001b[0m "
                  << fmt::format(fmt, std::forward<Args>(args)...) << std::endl;
    }

  private:
    LoggerSingleton() {}

    static LoggerSingleton &getInstanceImpl() {
        static LoggerSingleton instance{};
        return instance;
    }

    Level _level;
};
} // namespace detail

inline void SetVerbosityLevel(const Level &level) { detail::LoggerSingleton::getInstance().SetVerbosityLevel(level); }

template <typename... Args> void Info(fmt::format_string<Args...> fmt, Args &&...args) {
    detail::LoggerSingleton::getInstance().Info(fmt, std::forward<Args>(args)...);
}

template <typename... Args> void Warn(fmt::format_string<Args...> fmt, Args &&...args) {
    detail::LoggerSingleton::getInstance().Warn(fmt, std::forward<Args>(args)...);
}

template <typename... Args> void Error(fmt::format_string<Args...> fmt, Args &&...args) {
    detail::LoggerSingleton::getInstance().Error(fmt, std::forward<Args>(args)...);
}

template <typename... Args> void Debug(fmt::format_string<Args...> fmt, Args &&...args) {
    detail::LoggerSingleton::getInstance().Debug(fmt, std::forward<Args>(args)...);
}

} // namespace plasmatic::Log
