#pragma once

#include "Log.h"

namespace plasmatic {
template <typename... Args> [[noreturn]] void Abort(fmt::format_string<Args...> fmt, Args &&...args) {
    Log::Error(fmt, std::forward<Args>(args)...);

    std::abort();
}
} // namespace plasmatic
