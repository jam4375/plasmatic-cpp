#include "Mesh/Mesh.h"
#include "Utility/Utility.h"

#include <cxxopts.hpp>

namespace plasmatic {
static auto Run(const std::string &command, const std::string &mesh_filepath) -> int {
    if (command == "surface_mesh") {
        Mesh mesh(mesh_filepath);
        mesh.WriteSurfaceMesh("surface_mesh");
    } else {
        Abort("Unknown command: {}", command);
    }

    return 0;
}
} // namespace plasmatic

auto main(int argc, char **argv) -> int {
    std::unique_ptr<cxxopts::Options> options;
    try {
        options = std::make_unique<cxxopts::Options>(
            "Plasmatic", "Plasmatic: finite element modeling and simulation (Version: ${${PROJECT_NAME}_GIT_VERSION})");

        // clang-format off
        options->add_options()
            ("h,help", "Print usage message")
            ("v,verbosity", "Logging verbosity level (Debug, Info, Warn, or Error), default is 'Debug'", cxxopts::value<std::string>())
            ("m,mesh", "Mesh file (*.msh)", cxxopts::value<std::string>())
            ("c,command", "Command (surface_mesh, ...)", cxxopts::value<std::string>())
        ;
        // clang-format on

        auto result = options->parse(argc, argv);

        if (result.count("help") == 1) {
            std::cout << options->help() << std::endl;
            std::exit(0); // NOLINT(concurrency-mt-unsafe)
        }

        if (result.count("verbosity") == 1) {
            auto verbosity = result["verbosity"].as<std::string>();

            if (verbosity == "Debug") {
                plasmatic::Log::SetVerbosityLevel(plasmatic::Log::Level::Debug);
            } else if (verbosity == "Info") {
                plasmatic::Log::SetVerbosityLevel(plasmatic::Log::Level::Info);
            } else if (verbosity == "Warn") {
                plasmatic::Log::SetVerbosityLevel(plasmatic::Log::Level::Warn);
            } else if (verbosity == "Error") {
                plasmatic::Log::SetVerbosityLevel(plasmatic::Log::Level::Error);
            } else {
                throw std::runtime_error("Unknown verbosity level");
            }
        }

        std::string mesh_filepath;
        if (result.count("mesh") == 1) {
            mesh_filepath = result["mesh"].as<std::string>();
        }

        std::string command;
        if (result.count("command") == 1) {
            command = result["command"].as<std::string>();
        }

        // Entry point:
        return plasmatic::Run(command, mesh_filepath);
    } catch (const cxxopts::option_not_exists_exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cout << std::endl;
        std::cout << options->help() << std::endl;
        std::exit(1); // NOLINT(concurrency-mt-unsafe)
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cout << std::endl;
        std::cout << options->help() << std::endl;
        std::exit(1); // NOLINT(concurrency-mt-unsafe)
    }
}
