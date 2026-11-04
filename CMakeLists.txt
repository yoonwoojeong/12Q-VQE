cmake_minimum_required(VERSION 3.16)
project(qchem-vqe VERSION 0.2.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(CURL REQUIRED)

set(LIBCINT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/libcint" CACHE PATH "Path to libcint prefix")

file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
list(FILTER LIB_SOURCES EXCLUDE REGEX "src/main\\.c$")

add_library(qchem STATIC ${LIB_SOURCES})
target_include_directories(qchem PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${LIBCINT_DIR}/include"
)
target_compile_definitions(qchem PUBLIC
    QC_GIT_REV="cmake"
    QC_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
)
target_link_directories(qchem PUBLIC "${LIBCINT_DIR}/lib")
target_link_libraries(qchem PUBLIC
    OpenMP::OpenMP_C
    cint
    quadmath
    lapacke
    openblas
    nlopt
    cjson
    CURL::libcurl
    m
)

add_executable(qchem-vqe "${CMAKE_CURRENT_SOURCE_DIR}/src/main.c")
target_link_libraries(qchem-vqe PRIVATE qchem)
