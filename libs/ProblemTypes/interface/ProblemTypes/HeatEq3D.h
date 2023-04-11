#pragma once

#include "Mesh/Mesh.h"

#include <filesystem>

namespace plasmatic {

class HeatEq3D {
  public:
    struct Input {
        std::filesystem::path mesh_filename;
        Float thermal_conductivity = std::numeric_limits<Float>::quiet_NaN();
        std::unordered_map<std::string, Float> dirichlet_bcs = {};
        std::unordered_map<std::string, Float> neumann_bcs = {};
    };

    HeatEq3D(const Input &input);

    void Solve();

    void WriteVTK(const std::filesystem::path &output_filename) { _mesh.WriteVTK(output_filename); }

  private:
    Input _input;
    Mesh _mesh;
};

} // namespace plasmatic
