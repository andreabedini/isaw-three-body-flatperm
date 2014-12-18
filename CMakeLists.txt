cmake_minimum_required(VERSION 3.0)
project (isaw-three-body)

set(Boost_USE_STATIC_LIBS    on) 
set(Boost_USE_STATIC_RUNTIME on)
find_package(Boost REQUIRED COMPONENTS date_time program_options system thread)

set(HDF5_USE_STATIC_LIBRARIES on)
find_package(HDF5 REQUIRED COMPONENTS C HL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic -Wextra")

add_definitions(-DPACKAGE="${CMAKE_PROJECT_NAME}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_executable(main main.cpp)
target_link_libraries(main ${Boost_LIBRARIES} ${HDF5_LIBRARIES})

