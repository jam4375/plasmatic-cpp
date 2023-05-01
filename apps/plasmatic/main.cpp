#include "ProblemTypes/ProblemTypes.h"

#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <fstream>

namespace plasmatic {
static auto Run(const nlohmann::json &input) -> int {
    auto command = input["command"].get<std::string>();

    if (command == "surface_mesh") {
        auto mesh_filepath = input["mesh_filepath"].get<std::string>();
        Mesh mesh(mesh_filepath);
        mesh.WriteSurfaceMesh("surface_mesh");
    } else if (command == "run_thermal_sim") {
        HeatEq2D::Input thermal_input = {.mesh_filename = input["mesh_filepath"].get<std::string>(),
                                         .thermal_conductivity = input["thermal_conductivity"].get<Float>(),
                                         .dirichlet_bcs = {{}},
                                         .neumann_bcs = {{}}};

        for (const auto &item : input["dirichlet_bcs"].items()) {
            thermal_input.dirichlet_bcs.insert(
                {item.value()["surface_name"].get<std::string>(), item.value()["value"].get<Float>()});
        }

        for (const auto &item : input["neumann_bcs"].items()) {
            thermal_input.neumann_bcs.insert(
                {item.value()["surface_name"].get<std::string>(), item.value()["value"].get<Float>()});
        }

        HeatEq2D problem(thermal_input);

        problem.Solve();

        problem.WriteVTK(input["output_file"].get<std::string>() + ".vtk");
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
            ("i,input", "JSON input file", cxxopts::value<std::string>())
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

        nlohmann::json input = {};
        if (result.count("input") == 1) {
            std::ifstream in(result["input"].as<std::string>());
            input = nlohmann::json::parse(in);
            in.close();
        }

        // Entry point:
        return plasmatic::Run(input);
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
