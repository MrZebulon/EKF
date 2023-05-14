cmake_minimum_required(VERSION 3.22)
project(shai)

set(CMAKE_CXX_STANDARD 23)

file(GLOB_RECURSE SOURCES "src/*.cpp" "includes/*.hpp")
add_executable(shai ${SOURCES} src/models/TranslationRotationModel.cpp)
target_include_directories(shai PRIVATE includes SYSTEM eigen3)