#pragma once

#include "Element.h"
#include "Line.h"
#include "LineOrder2.h"
#include "Tetrahedron.h"
#include "TetrahedronOrder2.h"
#include "Triangle.h"
#include "TriangleOrder2.h"

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

    void AddTensorField(const std::string &field_name);

    std::shared_ptr<Element> GetElement(Integer dimension, Integer element_id) const {
        return _elements.at(static_cast<size_t>(dimension))[static_cast<size_t>(element_id)];
    }

    void ScalarFieldSetValue(const std::string &field_name, Integer index, Float value);

    void VectorFieldSetValue(const std::string &field_name, Integer index, std::array<Float, 3> value);

    void TensorFieldSetValue(const std::string &field_name, Integer index, std::array<Float, 6> value);

    const std::vector<Integer> &GetEntity(Integer dimension, Integer entity_id) const {
        return _entities.at(entity_id).at(static_cast<size_t>(dimension));
    }

    const std::vector<Integer> &GetPhysicalEntity(const std::string &physical_name, Integer dimension) const {
        return _physicalEntities.at(physical_name)[static_cast<size_t>(dimension)];
    }

    void WriteSurfaceMesh(const std::filesystem::path &base_filename) const;

  private:
    std::shared_ptr<std::vector<Coord>> _nodes;
    std::array<std::vector<std::shared_ptr<Element>>, 4> _elements;

    std::unordered_map<std::string, std::array<std::vector<Integer>, 4>> _physicalEntities;

    std::unordered_map<Integer, std::array<std::vector<Integer>, 4>> _entities;

    std::unordered_map<std::string, std::vector<Float>> _scalarFields;
    std::unordered_map<std::string, std::vector<std::array<Float, 3>>> _vectorFields;
    std::unordered_map<std::string, std::vector<std::array<Float, 6>>> _tensorFields;
};

} // namespace plasmatic
