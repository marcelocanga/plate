cmake_minimum_required(VERSION 3.8.2)

set(myapp plate)
project(${myapp}  CXX)
#
#*********************************************      system libraries
#

if(${U_ARCH} MATCHES  mac_64 )

 include_directories ( src/inc )

find_library(BLAS_LIB      NAMES blas      HINTS /usr/lib )
find_library(LAPACK_LIB    NAMES lapack    HINTS /usr/lib )

endif()

#
#*********************************************      package libraries
#

message("CMAKE_PREFIX_PATH:" ${CMAKE_PREFIX_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin) 
set(HEADER_FILES src/tool/Array.hh )
set(MOST src/main.cc src/Shape.cc src/Plate.cc src/Load.cc src/Support.cc src/Diagnostic.cc src/Report.cc src/Point.cc src/Solver.cc src/tool/Essential.cc )

add_executable(${myapp}  ${MOST} ${HEADER_FILES} )
target_include_directories(${myapp} PRIVATE src/tool)
target_link_libraries(${myapp}  -L${CMAKE_CURRENT_SOURCE_DIR}/lib/${U_ARCH} ${LAPACK_LIB} ${BLAS_LIB} )
set_target_properties(${myapp}  PROPERTIES LINKER_LANGUAGE CXX)

#
#*********************************************      
#
