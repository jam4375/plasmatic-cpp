#pragma once

#include "Coord.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Element {
  public:
    virtual ~Element();

    virtual Integer NumNodes() const = 0;

    virtual Integer VTKCellType() const = 0;

    virtual Float ShapeFn(Integer index, const Coord &coord) const = 0;

    virtual Integer GetNodeInex(Integer index) const = 0;

  private:
};

} // namespace plasmatic
