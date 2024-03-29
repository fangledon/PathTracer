cmake_minimum_required(VERSION 2.8)
project(pathtracer CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()


if (WIN32)
    set(embree_DIR "C:\\Program Files\\Intel\\Embree3 x64" CACHE PATH "Embree Directory")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")
    set(ARCH_LIBRARIES "")
else()
    if (APPLE)
        set(embree_DIR "/opt/local/lib/cmake/embree-3.8.0/" CACHE PATH "Embree Directory")
    elseif (UNIX)
        set(embree_DIR "/usr/lib64/cmake/embree-3.8.0/" CACHE PATH "Embree Directory")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -Wno-error=unused-variable -Wno-error=unused-function -Wno-error=unused-parameter")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    set(ARCH_LIBRARIES "pthread")
endif()

find_package(embree 3.0 REQUIRED)

set(LIB_SOURCES
    src/Engine.cpp
    src/Scene.cpp
    src/SceneLoader.cpp
    src/RenderPool.cpp
    src/RayTracerIntegrator.cpp
    src/lodepng/lodepng.cpp
)
add_definitions(-DGLM_ENABLE_EXPERIMENTAL)
include_directories(src/ ${EMBREE_INCLUDE_DIRS})
add_library(pathtracer_lib ${LIB_SOURCES})
set_target_properties(pathtracer_lib PROPERTIES OUTPUT_NAME pathtracer)
target_link_libraries(pathtracer_lib ${ARCH_LIBRARIES} ${EMBREE_LIBRARY})

set(CLI_SOURCES
    src/CliMain.cpp
)
add_executable(pathtracer_cli ${CLI_SOURCES})
set_target_properties(pathtracer_cli PROPERTIES OUTPUT_NAME pathtracer)
target_link_libraries(pathtracer_cli pathtracer_lib)

if (WIN32)
    set(WIN_SOURCES
        src/WindowsMain.cpp
    )
    add_executable(pathtracer_win ${WIN_SOURCES})
    target_link_libraries(pathtracer_win pathtracer_lib)
endif()
