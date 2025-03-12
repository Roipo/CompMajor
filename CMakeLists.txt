cmake_minimum_required(VERSION 3.0)
project(output)

# Option to enable PARDISO
option(USE_PARDISO "Use Intel MKL PARDISO solver" OFF)

# Set libigl paths
set(LIBIGL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libigl/include)
set(LIBIGL_EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libigl/external)

# Add Eigen include directory
include_directories(
    ${LIBIGL_INCLUDE_DIR}
    ${LIBIGL_EXTERNAL_DIR}/eigen
)

option(LIBIGL_WITH_OPENGL      "Use OpenGL"         OFF)

# Your source files
set(SRCFILES
    src/Energy.cpp
    src/EnergySymDir.cpp
    src/Newton.cpp
    src/Solver.cpp
    src/Utils.cpp
    src/main_cmd.cpp
    src/PardisoSolver.cpp
)


set(MKL_ARCH "intel64")
set(MKL_LINK "static")
set(MKL_THREADING "sequential")
find_package(MKL REQUIRED)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})



# Create executable
add_executable(${PROJECT_NAME}_bin ${SRCFILES})

# Set C++11 standard
set_target_properties(${PROJECT_NAME}_bin PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)


IF(USE_PARDISO)
    message("Use PARDISO")
    set(SRCFILES ${SRCFILES} src/PardisoSolver.cpp)
ENDIF()

target_link_libraries(${PROJECT_NAME}_bin PUBLIC MKL::MKL)
