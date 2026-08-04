#----------------------------------------------------------------------------
# A complete, runnable Geant4 example for UCMuGen.
#
#   cmake -S . -B build && cmake --build build -j
#   ./build/ucmugen_example run.mac        # batch
#   ./build/ucmugen_example                # interactive, with visualisation
#
# UCMuGen itself needs nothing built: it is a header. The only reason this
# directory has a build system at all is to give you something that runs.
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)
project(ucmugen_example)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

add_executable(ucmugen_example
  ucmugen_example.cc
  PrimaryGeneratorAction.cc
)

target_include_directories(ucmugen_example PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include      # UCMuGen.h
)

target_compile_features(ucmugen_example PRIVATE cxx_std_17)
target_link_libraries(ucmugen_example PRIVATE ${Geant4_LIBRARIES})

# Optional: the location- and date-aware PARMA flux. It is a separate header
# because it carries JAEA's non-commercial licence, so it is opt-in here too.
#   cmake -S . -B build -DUCMUGEN_PARMA=ON
option(UCMUGEN_PARMA "Use the PARMA/EXPACS flux (JAEA, non-commercial only)" OFF)
if(UCMUGEN_PARMA)
  target_compile_definitions(ucmugen_example PRIVATE UCMUGEN_EXAMPLE_PARMA)
  message(STATUS "PARMA flux enabled: this build is NON-COMMERCIAL ONLY")
endif()

# Copy the macros next to the binary so the example runs from the build dir.
configure_file(run.mac ${CMAKE_CURRENT_BINARY_DIR}/run.mac COPYONLY)
configure_file(vis.mac ${CMAKE_CURRENT_BINARY_DIR}/vis.mac COPYONLY)
configure_file(tracks.mac ${CMAKE_CURRENT_BINARY_DIR}/tracks.mac COPYONLY)
