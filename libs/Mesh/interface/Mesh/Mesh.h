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

    void WriteVTK(const std::filesystem::path &filename) const;

    Integer GetNumNodes() const { return static_cast<Integer>(_nodes->size()); }

    Integer GetNumElements(Integer dimension) const {
        return static_cast<Integer>(_elements.at(static_cast<size_t>(dimension)).size());
    }

    void AddScalarField(const std::string &field_name);

    void AddVectorField(const std::string &field_name);

    std::shared_ptr<Element> GetElement(Integer dimension, Integer element_id) {
        return _elements.at(static_cast<size_t>(dimension))[static_cast<size_t>(element_id)];
    }

  private:
    std::shared_ptr<std::vector<Coord>> _nodes;
    std::array<std::vector<std::shared_ptr<Element>>, 4> _elements;

    std::unordered_map<Integer, std::vector<Integer>> _entities;

    std::unordered_map<std::string, std::vector<Float>> _scalarFields;
    std::unordered_map<std::string, std::vector<std::array<Float, 3>>> _vectorFields;
};

} // namespace plasmatic
