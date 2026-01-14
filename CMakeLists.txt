# CMakeLists.txt
cmake_minimum_required(VERSION 4.0)  # Minimum CMake version

# Enable generation of compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

project(cppFFT)                      # Project name

# Set C++ standard (e.g., C++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # If c++17 isn't available, doesn't compile

# Add source files
set(SOURCES
    src/main.cpp
    src/fft.cpp
)

# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create an executable
add_executable(out ${SOURCES})
