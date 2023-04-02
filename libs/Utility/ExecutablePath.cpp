#include "interface/Utility/ExecutablePath.h"

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <sstream>
#include <unistd.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

#include <cstdint>
#include <vector>

namespace plasmatic {

std::filesystem::path GetExecutablePath() {
    uint32_t bufferSize = 512; // NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    std::vector<char> buffer(bufferSize + 1);

#if defined(_WIN32)
    ::GetModuleFileName(NULL, buffer.data(), bufferSize);

#elif defined(__linux__)
    // Get the process ID.
    int pid = getpid();

    // Construct a path to the symbolic link pointing to the process executable.
    // This is at /proc/<pid>/exe on Linux systems (we hope).
    std::ostringstream oss;
    oss << "/proc/" << pid << "/exe";
    std::string link = oss.str();

    // Read the contents of the link.
    int count = readlink(link.c_str(), buffer.data(), bufferSize);
    Check(count != -1, "Could not read symbolic link");
    buffer[count] = '\0';

#elif defined(__APPLE__)
    if (_NSGetExecutablePath(buffer.data(), &bufferSize) != 0) {
        buffer.resize(bufferSize);
        _NSGetExecutablePath(buffer.data(), &bufferSize);
    }

#else
#error Cannot yet find the executable on this platform
#endif

    std::string s = buffer.data();
    return std::filesystem::path(s).remove_filename();
}

} // namespace plasmatic
