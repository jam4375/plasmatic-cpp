#include "ExampleLib/ExampleLib.h"

#include <cxxopts.hpp>

namespace templo::ExampleApp {
static auto Run() -> int {
    ExampleLib::PrintMessage("Hello world");

    return 0;
}
} // namespace templo::ExampleApp

auto main(int argc, char **argv) -> int {
    std::unique_ptr<cxxopts::Options> options;
    try {
        options =
            std::make_unique<cxxopts::Options>("ExampleApp", "ExampleApp: example application for the templo template "
                                                             "repository (Version: ${${PROJECT_NAME}_GIT_VERSION})");

        // clang-format off
        options->add_options()
            ("h,help", "Print usage message")
        ;
        // clang-format on

        auto result = options->parse(argc, argv);

        if (result.count("help") == 1) {
            std::cout << options->help() << std::endl;
            std::exit(0); // NOLINT(concurrency-mt-unsafe)
        }

        // Entry point:
        return templo::ExampleApp::Run();
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