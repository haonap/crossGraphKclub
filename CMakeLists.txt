cmake_minimum_required(VERSION 3.16)
project(crossGraphKclub)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -std=c++11 -g")
# declara header gurobi directory
set(GUROBI_INCLUDE_DIRECTORY "/Library/gurobi951/macos_universal2/include")
# adicionar header gurobi
include_directories(${GUROBI_INCLUDE_DIRECTORY})
# declara library gurobi directory
set(GUROBI_LIB_DIRECTORY "/Library/gurobi951/macos_universal2/lib")
# adicionar library gurobi
link_directories(${GUROBI_LIB_DIRECTORY})
# criar executavel
set(SOURCE_FILES  main.cpp functions.cpp classes.cpp  ${GUROBI_INCLUDE_DIRECTORY}/gurobi_c++.h functions.h classes.h)
add_executable(crossGraphKclub ${SOURCE_FILES})

# adicionar libraries
target_link_libraries (crossGraphKclub gurobi_c++)
target_link_libraries (crossGraphKclub gurobi95)
