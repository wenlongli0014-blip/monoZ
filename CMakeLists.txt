cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# By default, set the build type to Release
if(NOT DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(default_build_type Release)
  message("Using default build type ${default_build_type}")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Select build type" FORCE
  )
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS "" Debug Release MinSizeRel RelWithDebInfo
  )
endif()

project(hzz2l2nu CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")

# Require C++17 in all targets
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -pedantic)

# Find external dependencies
find_package(
  Boost 1.72 REQUIRED
  COMPONENTS log program_options stacktrace_basic
)
find_package(ROOT 6 REQUIRED)
find_package(yaml-cpp REQUIRED)


# Minimalistic library to embed the hash of the current Git commit into the
# binary code as its version. CMake will be rerun whenever .git/index changes.
set_property(DIRECTORY APPEND
  PROPERTY CMAKE_CONFIGURE_DEPENDS
  .git/index
)

find_package(Git QUIET REQUIRED)
execute_process(
  COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(src/Version.cc.in "${CMAKE_BINARY_DIR}/Version.cc")

add_library(version STATIC
  "${CMAKE_BINARY_DIR}/Version.cc"
)
target_include_directories(version PRIVATE include)


# Jet calibration library
add_library(jerc STATIC
  src/JERC/FactorizedJetCorrector.cc
  src/JERC/JetCorrectionUncertainty.cc
  src/JERC/JetCorrectorParameters.cc
  src/JERC/JetResolution.cc
  src/JERC/JetResolutionObject.cc
  src/JERC/SimpleJetCorrectionUncertainty.cc
  src/JERC/SimpleJetCorrector.cc
)
target_include_directories(jerc PRIVATE src/JERC)
target_link_libraries(jerc
  PRIVATE ROOT::GenVector ROOT::Hist
)

# Library to access b tagging scale factors
add_library(btag STATIC
  src/BTag/BTagCalibrationStandalone.cc
)
target_include_directories(btag PRIVATE src/BTag)
target_link_libraries(btag
  PRIVATE ROOT::Hist
)

# Manually request that auxiliary libraries above are compiled with option
# POSITION_INDEPENDENT_CODE so that they can be included in the main shared
# library below. In newer versions of CMake could use object libraries instead.
set_target_properties(btag jerc version
  PROPERTIES POSITION_INDEPENDENT_CODE ON
)

add_library(hzz2l2nu SHARED
  src/AnalysisCommon.cc
  src/BTagger.cc
  src/BTagWeight.cc
  src/CollectionBuilder.cc
  src/Dataset.cc
  src/DileptonTrees.cc
  src/EGammaFromMisid.cc
  src/ElectronBuilder.cc
  src/ElectronTrees.cc
  src/EventTrees.cc
  src/EWCorrectionWeight.cc
  src/FileInPath.cc
  src/GenJetBuilder.cc
  src/GenPhotonBuilder.cc
  src/GenWeight.cc
  src/GenZZBuilder.cc
  src/GJetsWeight.cc
  src/JetBuilder.cc
  src/JetCorrector.cc
  src/JetGeometricVeto.cc
  src/KFactorCorrection.cc
  src/L1TPrefiringWeight.cc
  src/LeptonWeight.cc
  src/Logger.cc
  src/MeKinFilter.cc
  src/MetFilters.cc
  src/MetXYCorrections.cc
  src/MuonBuilder.cc
  src/NrbAnalysis.cc
  src/NrbTrees.cc
  src/Options.cc
  src/PhotonBuilder.cc
  src/PhotonPrescales.cc
  src/PhotonTrees.cc
  src/PhotonWeight.cc
  src/IsoTrackBuilder.cc
  src/PileUpIdFilter.cc
  src/PileUpIdWeight.cc
  src/PileUpWeight.cc
  src/PSWeight.cc
  src/PtMissBuilder.cc
  src/RoccoR.cc
  src/RunSampler.cc
  src/SmartSelectionMonitor.cc
  src/SmartSelectionMonitor_hzz.cc
  src/EventNumberFilter.cc
  src/TabulatedRandomGenerator.cc
  src/TauBuilder.cc
  src/TriggerFilter.cc
  src/TriggerWeight.cc
  src/TrileptonTrees.cc
  src/Utils.cc
  src/WeightCollector.cc
  src/ZGammaTrees.cc
)
target_include_directories(hzz2l2nu PUBLIC include)
target_link_libraries(hzz2l2nu
  PRIVATE stdc++fs
  PRIVATE jerc btag
  PRIVATE version
  PUBLIC Boost::boost Boost::log Boost::program_options
  PUBLIC Boost::stacktrace_basic -rdynamic  # To preserve human-readable names
  PUBLIC ROOT::Hist ROOT::MathCore ROOT::Physics ROOT::Tree ROOT::TreePlayer
  PUBLIC yaml-cpp
)

add_executable(runHZZanalysis src/runHZZanalysis.cc)
target_link_libraries(runHZZanalysis PRIVATE hzz2l2nu Boost::boost)
add_executable (nrbTreeHandler 
  src/nrbTreeHandler.cc 
  src/FileInPath.cc
  src/Logger.cc
  src/Options.cc
)
target_include_directories(nrbTreeHandler PUBLIC include)
target_link_libraries(nrbTreeHandler 
  PRIVATE stdc++fs
  PRIVATE version
  PUBLIC Boost::boost Boost::log Boost::program_options
  PUBLIC Boost::stacktrace_basic -rdynamic
  PUBLIC ROOT::Hist ROOT::ROOTDataFrame
  PUBLIC yaml-cpp
)
# Export targets for dependent projects
export(TARGETS hzz2l2nu
  NAMESPACE hzz2l2nu::
  FILE "${CMAKE_SOURCE_DIR}/cmake/hzz2l2nuTargets.cmake"
)

