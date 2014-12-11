project(exdepth)

cmake_minimum_required(VERSION 2.8)

find_package(Boost REQUIRED)

include_directories(${BOOST_INCLUDES})

add_executable(deep deep.cpp)
target_link_libraries(deep boost_program_options)

