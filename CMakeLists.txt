cmake_minimum_required(VERSION 3.12)
project(GrapeneTopology VERSION 0.1.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

set(CMAKE_USE_RELATIVE_PATHS ON)
set(CMAKE_SUPPRESS_REGENERATION ON)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(grapology)

target_include_directories(grapology PRIVATE ./ ./include)

file(GLOB SOURCES "src/*.cpp" "src/*.cxx" "src/*.cc")
target_sources(grapology PRIVATE ${SOURCES})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(grapology PRIVATE
        $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>
        $<$<CONFIG:MinSizeRel>:-Os -march=native -DNDEBUG>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g -march=native>
        $<$<CONFIG:Debug>:-O0 -g -DDEBUG>
    )
elseif(MSVC)
    target_compile_options(grapology PRIVATE
        $<$<CONFIG:Release>:/O2 /DNDEBUG>
        $<$<CONFIG:MinSizeRel>:/O1 /DNDEBUG>
        $<$<CONFIG:RelWithDebInfo>:/O2 /Zi>
        $<$<CONFIG:Debug>:/Od /Zi /DDEBUG>
    )
else()
    target_compile_options(grapology PRIVATE
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
        $<$<CONFIG:Debug>:-O0 -g -DDEBUG>
    )
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Using ${CMAKE_CXX_STANDARD} standard")