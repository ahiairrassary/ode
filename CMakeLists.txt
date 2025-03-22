cmake_minimum_required(VERSION 3.22)

project(ODE LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(USER_COMPILER_WARNINGS -Wall -Wextra -Werror -Wfatal-errors -pedantic -pedantic-errors -Wshadow)

############################################################
## ODE
############################################################

set(SRC_FILES
    src/main.cpp
)

set(INCLUDE_FILES
    src/OdeSolver.hpp
)

add_executable(ODE
    ${SRC_FILES}
    ${INCLUDE_FILES}
)

target_compile_options(ODE
    PRIVATE
        ${USER_COMPILER_WARNINGS}
)

target_include_directories(ODE
    PRIVATE
        src
)
