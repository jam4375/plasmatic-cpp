#pragma once

#include "Utility/Utility.h"

namespace plasmatic {

struct Coord {
    Float x = std::numeric_limits<Float>::quiet_NaN();
    Float y = std::numeric_limits<Float>::quiet_NaN();
    Float z = std::numeric_limits<Float>::quiet_NaN();
};

} // namespace plasmatic
