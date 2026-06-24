cmake_minimum_required(VERSION 3.16)
project(MuonRock)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

option(USE_ECOMUG "Compile with EcoMug generator" OFF)
if(USE_ECOMUG)
    # If ECOMUG_INCLUDE is provided, add it; otherwise assume EcoMug.h
    # is already in the project's own include/ directory.
    if(DEFINED ECOMUG_INCLUDE)
        include_directories(${ECOMUG_INCLUDE})
        message(STATUS "EcoMug enabled — include dir: ${ECOMUG_INCLUDE}")
    else()
        message(STATUS "EcoMug enabled — using project include/ directory")
    endif()
    add_compile_definitions(USE_ECOMUG)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(MuonRock MuonRock.cc ${sources} ${headers})
target_link_libraries(MuonRock ${Geant4_LIBRARIES})

foreach(_mac run_file.mac run_powerlaw.mac run_ecomug.mac vis.mac)
    if(EXISTS "${PROJECT_SOURCE_DIR}/${_mac}")
        configure_file(${PROJECT_SOURCE_DIR}/${_mac}
                       ${PROJECT_BINARY_DIR}/${_mac} COPYONLY)
    endif()
endforeach()
