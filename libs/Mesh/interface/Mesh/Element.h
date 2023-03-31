#pragma once

#include "Coord.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Element {
  public:
    virtual ~Element();

    virtual Float ShapeFn(Integer index, const Coord &coord) = 0;

  private:
};

} // namespace plasmatic
