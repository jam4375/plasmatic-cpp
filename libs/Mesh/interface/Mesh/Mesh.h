#pragma once

#include "Element.h"
#include "Triangle.h"

#include "Utility/Utility.h"

#include <array>
#include <filesystem>
#include <vector>

namespace plasmatic {

class Mesh {
  public:
    Mesh(const std::filesystem::path &filename);

  private:
    std::shared_ptr<std::vector<Coord>> _nodes;
    std::vector<std::shared_ptr<Element>> _elements;
};

} // namespace plasmatic
