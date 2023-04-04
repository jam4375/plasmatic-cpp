#pragma once

#include "Mesh/Mesh.h"

#include <filesystem>

namespace plasmatic {

class ProblemType {
  public:
    ProblemType(const std::filesystem::path &mesh_filename);

    void Solve();

    void WriteVTK(const std::filesystem::path &output_filename) { _mesh.WriteVTK(output_filename); }

  private:
    Mesh _mesh;
};

} // namespace plasmatic
