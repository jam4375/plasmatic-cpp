#pragma once

#include "Element.h"
#include "Line.h"
#include "Triangle.h"

#include "Utility/Utility.h"

#include <array>
#include <filesystem>
#include <unordered_map>
#include <vector>

namespace plasmatic {

class Mesh {
  public:
    Mesh(const std::filesystem::path &filename);

    Integer GetNumNodes() const { return static_cast<Integer>(_nodes->size()); }

    Integer GetNumElements(Integer dimension) const {
        return static_cast<Integer>(_elements.at(static_cast<size_t>(dimension)).size());
    }

  private:
    std::shared_ptr<std::vector<Coord>> _nodes;
    std::array<std::vector<std::shared_ptr<Element>>, 4> _elements;

    std::unordered_map<Integer, std::vector<Integer>> _entities;
};

} // namespace plasmatic
