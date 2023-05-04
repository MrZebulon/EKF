cmake_minimum_required(VERSION 3.22)
project(shai)

set(CMAKE_CXX_STANDARD 23)

add_executable(shai src/main.cpp
        includes/BaseModel.hpp
        includes/TranslationRotationModel.hpp src/TranslationRotationModel.cpp
        includes/MEKF.hpp includes/SensorParameters.hpp includes/Utils.hpp)
