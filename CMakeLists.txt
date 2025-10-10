cmake_minimum_required(VERSION 2.8)
project(LICO)
#add_compile_options(-mavx2)
add_compile_options(-mlzcnt)
#add_compile_options(-msse4.2)
add_compile_options(-mbmi2)
add_compile_options(-mavx512f)
add_compile_options(-mavx512vl)  # Add this line
add_compile_options(-mavx512cd) # Add this line
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # Disable all warnings

  if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE "Release")
  endif()
  MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

set(CMAKE_CXX_FLAGS_NDEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")


set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external/Vp2SetOperation/include)

# sdsl-lite
include_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/include)
link_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/lib)

# FastPFor
include_directories(${PROJECT_SOURCE_DIR}/external/FastPFor/headers)
link_directories(${PROJECT_SOURCE_DIR}/external/FastPFor/build)

function(add_variants target_base source_file)
  set(target_z "${target_base}")
  add_executable(${target_z} ${source_file})
  set_target_properties(${target_z} PROPERTIES
    COMPILE_DEFINITIONS "RESIDUAL_COMPRESS=1")
  target_link_libraries(${target_z} sdsl FastPFOR)

  set(target_p "${target_base}_none")
  add_executable(${target_p} ${source_file})
  set_target_properties(${target_p} PROPERTIES
    COMPILE_DEFINITIONS "RESIDUAL_COMPRESS=0")
  target_link_libraries(${target_p} sdsl)
endfunction()

function(add_uncompress target_base source_file)
  set(target_p "${target_base}")
  add_executable(${target_p} ${source_file})
  set_target_properties(${target_p} PROPERTIES
          COMPILE_DEFINITIONS "RESIDUAL_COMPRESS=0")
  target_link_libraries(${target_p} sdsl FastPFOR)
endfunction()




add_variants(lico_build lico_build.cpp)

add_variants(lico_decode lico_decode.cpp)

add_variants(lico_query lico_query.cpp)

#add_executable(lico_estimate_none lico_estimate.cpp)
#set_target_properties(lico_estimate_none PROPERTIES
#        COMPILE_DEFINITIONS "RESIDUAL_COMPRESS=0")
#target_link_libraries(lico_estimate_none sdsl)
add_uncompress(lico_estimate lico_estimate.cpp)

add_uncompress(other other_index.cpp)
#add_executable(other other_index.cpp)
#target_link_libraries(other sdsl FastPFOR)