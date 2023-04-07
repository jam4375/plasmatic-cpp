#include "interface/Mesh/Mesh.h"

#include <fstream>
#include <sstream>

namespace plasmatic {

Mesh::Mesh(const std::filesystem::path &filename) : _nodes(std::make_shared<std::vector<Coord>>()) {
    Log::Info("Reading mesh from file '{}'", filename.string());
    std::ifstream in(filename);

    std::string line;

    // Maps physical tags to physical names
    std::unordered_map<Integer, std::string> physical_tags;

    while (std::getline(in, line)) {

        if (line == "$PhysicalNames") {
            Integer num_physical_entities = 0;

            {
                std::getline(in, line);
                std::stringstream ss(line);
                ss >> num_physical_entities;
            }

            for (Integer ii = 0; ii < num_physical_entities; ++ii) {
                Integer physical_dim = 0;
                Integer physical_tag = 0;
                std::string name;
                {
                    std::getline(in, line);
                    std::stringstream ss(line);
                    ss >> physical_dim >> physical_tag >> name;
                }

                // Remove leading and trailing quotes:
                name.erase(0, 1);
                name.erase(name.size() - 1);

                physical_tags.insert({physical_tag, name});
            }
        }

        if (line == "$Entities") {
            Integer num_points = 0;
            Integer num_curves = 0;
            Integer num_surfaces = 0;
            Integer num_volumes = 0;

            {
                std::getline(in, line);
                std::stringstream ss(line);
                ss >> num_points >> num_curves >> num_surfaces >> num_volumes;
            }

            for (Integer ii = 0; ii < num_points; ++ii) {
                Integer point_tag = 0;
                Float x = 0;
                Float y = 0;
                Float z = 0;
                Integer num_physical_tags = 0;

                std::getline(in, line);
                std::stringstream ss(line);
                ss >> point_tag >> x >> y >> z >> num_physical_tags;

                for (Integer jj = 0; jj < num_physical_tags; ++jj) {
                    Integer physical_tag = 0;
                    ss >> physical_tag;

                    _physicalEntities[physical_tags.at(physical_tag)][0].push_back(point_tag);
                }
            }

            for (Integer ii = 0; ii < num_curves; ++ii) {
                Integer curve_tag = 0;
                Float min_x = 0;
                Float min_y = 0;
                Float min_z = 0;
                Float max_x = 0;
                Float max_y = 0;
                Float max_z = 0;
                Integer num_physical_tags = 0;

                std::getline(in, line);
                std::stringstream ss(line);
                ss >> curve_tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z >> num_physical_tags;

                for (Integer jj = 0; jj < num_physical_tags; ++jj) {
                    Integer physical_tag = 0;
                    ss >> physical_tag;

                    _physicalEntities[physical_tags.at(physical_tag)][1].push_back(curve_tag);
                }
            }

            for (Integer ii = 0; ii < num_surfaces; ++ii) {
                Integer surface_tag = 0;
                Float min_x = 0;
                Float min_y = 0;
                Float min_z = 0;
                Float max_x = 0;
                Float max_y = 0;
                Float max_z = 0;
                Integer num_physical_tags = 0;

                std::getline(in, line);
                std::stringstream ss(line);
                ss >> surface_tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z >> num_physical_tags;

                for (Integer jj = 0; jj < num_physical_tags; ++jj) {
                    Integer physical_tag = 0;
                    ss >> physical_tag;

                    _physicalEntities[physical_tags.at(physical_tag)][2].push_back(surface_tag);
                }
            }

            for (Integer ii = 0; ii < num_volumes; ++ii) {
                Integer volume_tag = 0;
                Float min_x = 0;
                Float min_y = 0;
                Float min_z = 0;
                Float max_x = 0;
                Float max_y = 0;
                Float max_z = 0;
                Integer num_physical_tags = 0;

                std::getline(in, line);
                std::stringstream ss(line);
                ss >> volume_tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z >> num_physical_tags;

                for (Integer jj = 0; jj < num_physical_tags; ++jj) {
                    Integer physical_tag = 0;
                    ss >> physical_tag;

                    _physicalEntities[physical_tags.at(physical_tag)][3].push_back(volume_tag);
                }
            }
        }

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

                    constexpr auto dimension = 0; // 0d

                    _entities[entity_tag][static_cast<size_t>(dimension)].push_back(
                        static_cast<Integer>(_nodes->size()));
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

                    if (element_type == 1) {
                        // 2-node line

                        Integer element_tag = 0;
                        Integer ind1 = 0;
                        Integer ind2 = 0;

                        std::stringstream ss(line);
                        ss >> element_tag >> ind1 >> ind2;

                        std::array<Integer, 2> node_indices = {ind1 - 1, ind2 - 1};

                        constexpr auto dimension = 1; // 1d
                        _entities[entity_tag][static_cast<size_t>(dimension)].push_back(
                            static_cast<Integer>(_elements[static_cast<size_t>(dimension)].size()));
                        _elements[static_cast<size_t>(dimension)].push_back(
                            std::make_shared<Line>(node_indices, _nodes));
                    } else if (element_type == 2) {
                        // 3-node triangle

                        Integer element_tag = 0;
                        Integer ind1 = 0;
                        Integer ind2 = 0;
                        Integer ind3 = 0;

                        std::stringstream ss(line);
                        ss >> element_tag >> ind1 >> ind2 >> ind3;

                        std::array<Integer, 3> node_indices = {ind1 - 1, ind2 - 1, ind3 - 1};

                        constexpr auto dimension = 2; // 2d
                        _entities[entity_tag][static_cast<size_t>(dimension)].push_back(
                            static_cast<Integer>(_elements[static_cast<size_t>(dimension)].size()));
                        _elements[static_cast<size_t>(dimension)].push_back(
                            std::make_shared<Triangle>(node_indices, _nodes));
                    } else if (element_type == 4) {
                        // 4-node tetrahedron

                        Integer element_tag = 0;
                        Integer ind1 = 0;
                        Integer ind2 = 0;
                        Integer ind3 = 0;
                        Integer ind4 = 0;

                        std::stringstream ss(line);
                        ss >> element_tag >> ind1 >> ind2 >> ind3 >> ind4;

                        std::array<Integer, 4> node_indices = {ind1 - 1, ind2 - 1, ind3 - 1, ind4 - 1};

                        constexpr auto dimension = 3; // 3d
                        _entities[entity_tag][static_cast<size_t>(dimension)].push_back(
                            static_cast<Integer>(_elements[static_cast<size_t>(dimension)].size()));
                        _elements[static_cast<size_t>(dimension)].push_back(
                            std::make_shared<Tetrahedron>(node_indices, _nodes));
                    }
                }
            }
        }
    }

    Log::Info("Num nodes read = {}", _nodes->size());
    for (size_t ii = 0; ii < 4; ++ii) {
        Log::Info("Num elements read = {} (dimension {})", _elements.at(ii).size(), ii);
    }

    in.close();
}

void Mesh::WriteVTK(const std::filesystem::path &filename) const {
    std::ofstream out(filename);

    constexpr auto float_precision = 16;

    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "Generated by Plasmatic" << std::endl;
    out << "ASCII" << std::endl;

    out << "DATASET UNSTRUCTURED_GRID" << std::endl;
    out << "POINTS " << _nodes->size() << " double" << std::endl;
    for (size_t ii = 0; ii < _nodes->size(); ++ii) {
        out << std::setprecision(float_precision) << (*_nodes)[ii].x << " ";
        out << std::setprecision(float_precision) << (*_nodes)[ii].y << " ";
        out << std::setprecision(float_precision) << (*_nodes)[ii].z << std::endl;
    }
    out << std::endl;

    Integer size_of_elements = 0;
    Integer num_elements = 0;
    for (const auto &elements : _elements) {
        for (const auto &element : elements) {
            size_of_elements += element->NumNodes() + 1;
            num_elements++;
        }
    }

    out << "CELLS " << num_elements << " " << size_of_elements << std::endl;
    for (const auto &elements : _elements) {
        for (const auto &element : elements) {
            out << element->NumNodes();
            for (Integer jj = 0; jj < element->NumNodes(); ++jj) {
                out << " " << element->GetNodeIndex(jj);
            }
            out << std::endl;
        }
    }
    out << std::endl;

    out << "CELL_TYPES " << num_elements << std::endl;
    for (const auto &elements : _elements) {
        for (const auto &element : elements) {
            out << element->VTKCellType() << std::endl;
        }
    }
    out << std::endl;

    out << "POINT_DATA " << _nodes->size() << std::endl;
    for (const auto &[data_name, values] : _scalarFields) {
        out << "SCALARS " << data_name << " double" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;
        for (const auto &value : values) {
            out << std::setprecision(float_precision) << value << std::endl;
        }
        out << std::endl;
    }

    for (const auto &[data_name, values] : _vectorFields) {
        out << "VECTORS " << data_name << " double" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;
        for (const auto &value : values) {
            out << std::setprecision(float_precision) << value[0] << " " << value[1] << " " << value[2] << std::endl;
        }
        out << std::endl;
    }

    out.close();
}

void Mesh::AddScalarField(const std::string &field_name) {
    _scalarFields.insert({field_name, std::vector<Float>(_nodes->size())});
}

void Mesh::AddVectorField(const std::string &field_name) {
    _vectorFields.insert({field_name, std::vector<std::array<Float, 3>>(_nodes->size())});
}

void Mesh::ScalarFieldSetValue(const std::string &field_name, Integer index, Float value) {
    _scalarFields.at(field_name)[static_cast<size_t>(index)] = value;
}

void Mesh::VectorFieldSetValue(const std::string &field_name, Integer index, std::array<Float, 3> value) {
    _vectorFields.at(field_name)[static_cast<size_t>(index)] = value;
}

} // namespace plasmatic
