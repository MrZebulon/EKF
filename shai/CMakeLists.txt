cmake_minimum_required(VERSION 3.22)
project(shai)

set(CMAKE_CXX_STANDARD 23)

add_executable(shai src/main.cpp
        includes/BaseModel.hpp
        includes/TranslationRotationModel.hpp src/TranslationRotationModel.cpp
        includes/EKFEngine.hpp includes/SensorParameters.hpp includes/Utils.hpp includes/TranslationRotationModelV2.hpp includes/BaseDataCollector.hpp)

target_include_directories(shai PRIVATE includes)
target_include_directories(shai PRIVATE eigen3)