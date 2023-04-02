#include "interface/Mesh/Mesh.h"

#include <fstream>
#include <sstream>

namespace plasmatic {

Mesh::Mesh(const std::filesystem::path &filename) : _nodes(std::make_shared<std::vector<Coord>>()) {
    Log::Info("Reading mesh from file '{}'", filename.string());
    std::ifstream in(filename);

    std::string line;

    while (std::getline(in, line)) {
        if (line == "$Nodes") {
            Integer num_entity_blocks = 0;
            Integer num_nodes = 0;

            {
                std::getline(in, line);
                std::stringstream ss(line);
                ss >> num_entity_blocks >> num_nodes;
            }

            for (Integer ii = 0; ii < num_entity_blocks; ++ii) {
                Integer entity_dim = 0;
                Integer entity_tag = 0;
                Integer parametric = 0;
                Integer nodes_in_block = 0;
                {
                    std::getline(in, line);
                    std::stringstream ss(line);
                    ss >> entity_dim >> entity_tag >> parametric >> nodes_in_block;
                }

                for (Integer jj = 0; jj < nodes_in_block; ++jj) {
                    std::getline(in, line); // node tags
                }

                for (Integer jj = 0; jj < nodes_in_block; ++jj) {
                    Coord coord;

                    std::getline(in, line); // node coords
                    std::stringstream ss(line);
                    ss >> coord.x >> coord.y >> coord.z;

                    _entities[entity_tag].push_back(static_cast<Integer>(_nodes->size()));
                    _nodes->push_back(coord);
                }
            }
        }

        if (line == "$Elements") {
            Integer num_entity_blocks = 0;
            Integer num_elements = 0;

            {
                std::getline(in, line);
                std::stringstream ss(line);
                ss >> num_entity_blocks >> num_elements;
            }

            for (Integer ii = 0; ii < num_entity_blocks; ++ii) {
                Integer entity_dim = 0;
                Integer entity_tag = 0;
                Integer element_type = 0;
                Integer elements_in_block = 0;
                {
                    std::getline(in, line);
                    std::stringstream ss(line);
                    ss >> entity_dim >> entity_tag >> element_type >> elements_in_block;
                }

                for (Integer jj = 0; jj < elements_in_block; ++jj) {
                    std::getline(in, line);

                    if (element_type == 2) {
                        // 3-node triangle

                        Integer element_tag = 0;
                        std::array<Integer, 3> node_indices = {};

                        std::stringstream ss(line);
                        ss >> element_tag >> node_indices[0] >> node_indices[1] >> node_indices[2];

                        _elements.push_back(std::make_shared<Triangle>(node_indices, _nodes));
                    }
                }
            }
        }
    }

    Log::Info("Num nodes read = {}", _nodes->size());
    Log::Info("Num elements read = {}", _elements.size());

    in.close();
}

} // namespace plasmatic
