#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Mesh {
  public:
  private:
    std::shared_ptr<std::vector<Coord>> _nodes;
    std::vector<Element> _elements;
};

} // namespace plasmatic
