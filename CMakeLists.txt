cmake_minimum_required(VERSION 2.8.9)
project(bruceRL)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_STANDARD 11)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")

#________________________________________________________________###
find_package(CGAL)
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  ADD_DEFINITIONS( "-DDEBUG_TRACE" )
  find_package(Eigen3 3.1.0)

    if (EIGEN3_FOUND)
        include( ${EIGEN3_USE_FILE} )
      endif()

endif()
#________________________________________________________________###

add_library(bruceRL SHARED ${SOURCES})
install(TARGETS bruceRL DESTINATION /usr/lib)



