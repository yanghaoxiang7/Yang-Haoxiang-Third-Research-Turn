cmake_minimum_required(VERSION 2.8.7)

#initialize project
set(proj_name a1_mass_spring)
project(${proj_name})

#set paths
if(NOT proj_path)
	set(proj_src_path ${PROJECT_SOURCE_DIR})
endif(NOT proj_path)

if(NOT root_path)
	set(root_path ${proj_src_path}/../..)
	add_definitions(-DROOT_PATH=\"${root_path}\")
endif(NOT root_path)

#include Eigen
include_directories(${root_path}/ext/eigen)

#include hearder and source files in /src
file(GLOB_RECURSE src_cpp ${root_path}/src/*.cpp)
file(GLOB_RECURSE src_h ${root_path}/src/*.h)
list(APPEND src_files ${src_cpp} ${src_h})
source_group("src" FILES ${src_cpp} ${src_h})
include_directories(${root_path}/src)
	
#include hearder and source files in /proj/a*
file(GLOB_RECURSE proj_cpp ${proj_src_path}/*.cpp)
file(GLOB_RECURSE proj_h ${proj_src_path}/*.h)
list(APPEND src_files ${proj_cpp} ${proj_h})
source_group("proj" FILES ${proj_cpp} ${proj_h})
include_directories(${proj_path}/src)

#include glm
include_directories(${root_path}/viewer/ext/glm)

#include freeglut and glew
if(WIN32)
	set(freeglut_src_path ${root_path}/viewer/ext/freeglut/include)
	set(freeglut_lib_path ${root_path}/viewer/ext/freeglut/lib/x64)
	include_directories(${freeglut_src_path})

	set(glew_libs ${freeglut_lib_path}/glew32.lib)
	list(APPEND lib_files ${glew_libs})
	set(glut_libs debug ${freeglut_lib_path}/freeglutd.lib optimized ${freeglut_lib_path}/freeglut.lib)
	message(STATUS ${glut_libs})
	list(APPEND lib_files ${glut_libs})
		
elseif(APPLE)  #freeglut and glew are installed on macos by "brew install freeglut glew"
	find_package(GLEW REQUIRED)
	if (GLEW_FOUND)
    		include_directories(${GLEW_INCLUDE_DIRS})
    		link_libraries(${GLEW_LIBRARIES})
	endif()

	find_package(GLUT REQUIRED)
	if (GLUT_FOUND)
    		include_directories(${GLUT_INCLUDE_DIRS})
    		link_libraries(${GLUT_LIBRARIES})
	endif()

	find_package(OpenGL REQUIRED)
    	include_directories(${OpenGL_INCLUDE_DIRS})
    	link_directories(${OpenGL_LIBRARY_DIRS})
    	add_definitions(${OpenGL_DEFINITIONS})

	list(APPEND lib_files ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

elseif(UNIX) #freeglut and glew are installed on linux by "sudo apt-get install freeglut3-dev libglew-dev"
	set(GCC_COVERAGE_COMPILE_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS} -lGL -lglut -lGLU -lGLEW")
	set(GCC_COVERAGE_LINK_FLAGS "${GCC_COVERAGE_LINK_FLAGS} -lGL -lglut -lGLU -lGLEW")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
	set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif(WIN32)

#include hearder and source files in /viewer/src
file(GLOB_RECURSE viewer_src_cpp ${root_path}/viewer/src/*.cpp)
file(GLOB_RECURSE viewer_src_h ${root_path}/viewer/src/*.h)
list(APPEND src_files ${viewer_src_cpp} ${viewer_src_h})
source_group("viewer-src" FILES ${viewer_src_cpp} ${viewer_src_h})
include_directories(${root_path}/viewer/src)

#set compiling flags
set(CMAKE_CXX_STANDARD 11)	#c++11
if(UNIX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")	#c++11
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")	#turn off sign-compare warning
endif(UNIX)
if(WIN32)
	add_definitions(-D_DISABLE_EXTENDED_ALIGNED_STORAGE)	#fix compiling issue for VS2017
endif(WIN32)

#add executable
add_executable(${proj_name} ${src_files})
target_link_libraries(${proj_name} ${lib_files})
