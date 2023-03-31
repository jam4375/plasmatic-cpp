#pragma once

#include "Element.h"
#include "Triangle.h"

#include "Utility/Utility.h"

#include <array>
#include <filesystem>
#include <vector>
#include <unordered_map>

namespace plasmatic {

class Mesh {
  public:
    Mesh(const std::filesystem::path &filename);

  private:
    std::shared_ptr<std::vector<Coord>> _nodes;
    std::vector<std::shared_ptr<Element>> _elements;

    std::unordered_map<Integer, std::vector<Integer>> _entities;
};

} // namespace plasmatic
