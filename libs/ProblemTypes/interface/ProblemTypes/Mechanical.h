#pragma once

#include "Mesh/Mesh.h"

#include <filesystem>

namespace plasmatic {

class Mechanical {
  public:
    struct Input {
        std::filesystem::path mesh_filename;
        Float youngs_modulus = std::numeric_limits<Float>::quiet_NaN();
        Float poisson_ratio = std::numeric_limits<Float>::quiet_NaN();
        std::unordered_map<std::string, std::array<Float, 3>> dirichlet_bcs = {};
        std::unordered_map<std::string, std::array<Float, 3>> neumann_bcs = {};
    };

    Mechanical(const Input &input);

    void Solve();

    void WriteVTK(const std::filesystem::path &output_filename) { _mesh.WriteVTK(output_filename); }

  private:
    Input _input;
    Mesh _mesh;
};

} // namespace plasmatic
