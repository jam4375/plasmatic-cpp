#pragma once

#include "Log.h"

namespace plasmatic {
template <typename... Args> void Check(bool condition, fmt::format_string<Args...> fmt, Args &&...args) {
    if (condition) {
        return;
    }

    Log::Error(fmt, std::forward<Args>(args)...);

    std::abort();
}
} // namespace plasmatic
