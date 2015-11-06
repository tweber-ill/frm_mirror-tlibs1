project(tlibs)
cmake_minimum_required(VERSION 3.0)

#set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_BUILD_TYPE Release)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ".")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

find_package(Boost REQUIRED COMPONENTS iostreams system filesystem regex python REQUIRED)
find_package(Threads REQUIRED)
find_package(Qt4 REQUIRED)
find_package(Freetype REQUIRED)
set(CMAKE_AUTOMOC OFF)
set(CMAKE_AUTOUIC OFF)


#message("Project source dir: ${PROJECT_SOURCE_DIR}")
#message("CMake source dir: ${CMAKE_SOURCE_DIR}")
#message("Binary dir: ${CMAKE_BINARY_DIR}")
message("Install path prefix: ${CMAKE_INSTALL_PREFIX}")


add_definitions(-DUSE_IOSTR -DUSE_BOOST_REX)
add_definitions(-DNO_LAPACK)
#add_definitions(-DNO_QT)
#add_definitions(-DUSE_FFTW)
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(${Boost_CXX_FLAGS})


include_directories("${PROJECT_SOURCE_DIR}" ".")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${QT_INCLUDE_DIR}")
include_directories("${FREETYPE_INCLUDE_DIRS}")


add_definitions(-std=c++11)
#add_definitions(-march=native -s)


add_library(tlibs SHARED
	file/file.cpp  file/loadinstr.cpp  file/loadtxt.cpp  file/recent.cpp  file/x3d.cpp  file/xml.cpp
	helper/log.cpp helper/flags.cpp
	math/fourier.cpp  math/linalg2.cpp  math/rand.cpp  # math/minuit.cpp  math/fftw.cpp
	net/tcp.cpp
	string/spec_char.cpp
	gfx/gl.cpp  gfx/gnuplot.cpp
	version.cpp
)


target_link_libraries(tlibs ${Boost_LIBRARIES} Threads::Threads
	Qt4::QtCore Qt4::QtGui
	${FREETYPE_LIBRARY}
)


install(TARGETS tlibs  DESTINATION lib)

install(FILES ${PROJECT_SOURCE_DIR}/*.h  DESTINATION include/tlibs)
install(FILES ${PROJECT_SOURCE_DIR}/file/*.h  DESTINATION include/tlibs/file)
install(FILES ${PROJECT_SOURCE_DIR}/helper/*.h  DESTINATION include/tlibs/helper)
install(FILES ${PROJECT_SOURCE_DIR}/math/*.h ${PROJECT_SOURCE_DIR}/math/*.hpp  DESTINATION include/tlibs/math)
install(FILES ${PROJECT_SOURCE_DIR}/net/*.h  DESTINATION include/tlibs/net)
install(FILES ${PROJECT_SOURCE_DIR}/string/*.h  DESTINATION include/tlibs/string)
install(FILES ${PROJECT_SOURCE_DIR}/gfx/*.h  DESTINATION include/tlibs/gfx)
install(FILES ${PROJECT_SOURCE_DIR}/time/*.h  DESTINATION include/tlibs/time)
install(FILES ${PROJECT_SOURCE_DIR}/cl/*.h  DESTINATION include/tlibs/cl)

install(FILES ${PROJECT_SOURCE_DIR}/COPYING ${PROJECT_SOURCE_DIR}/AUTHORS ${PROJECT_SOURCE_DIR}/LICENSES
        DESTINATION share/tlibs
)
