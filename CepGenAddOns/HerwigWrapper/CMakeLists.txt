#--- searching for ThePEG 2
set(THEPEG_DIRS $ENV{THEPEG_DIR} ${THEPEG_DIR} /usr /usr/local)
find_library(THEPEG ThePEG HINTS ${THEPEG_DIRS} PATH_SUFFIXES lib lib/ThePEG)
find_library(THEPEG_LHE LesHouches.so HINTS ${THEPEG_DIRS} PATH_SUFFIXES lib lib/ThePEG)
find_path(THEPEG_INCLUDE ThePEG HINTS ${THEPEG_DIRS} PATH_SUFFIXES include)
if(NOT THEPEG)
  return()
endif()

#--- searching for Herwig 7
set(HERWIG_DIRS $ENV{HERWIG7_DIR} ${HERWIG7_DIR} /usr /usr/local)
find_library(HERWIG Herwig.so HINTS ${HERWIG_DIRS} PATH_SUFFIXES lib/Herwig)
find_library(HERWIG_API HerwigAPI HINTS ${HERWIG_DIRS} PATH_SUFFIXES lib/Herwig)
find_path(HERWIG_INCLUDE Herwig HINTS ${HERWIG_DIRS} PATH_SUFFIXES include include/Herwig)
if(NOT HERWIG OR NOT HERWIG_API)
  return()
endif()

#--- searching for the version of HepMC used for Herwig compilation
set(HEPMC_DIRS $ENV{HERWIG_HEPMC_DIR} /usr /usr/local)
find_library(HEPMC_LIB NAME HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES lib64 lib)
find_path(HEPMC_INCLUDE NAME HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES include)
if(NOT HEPMC_LIB)
  message(STATUS "Herwig was found, but the version of HepMC is incompatible.")
  return()
endif()
message(STATUS "Herwig API found in ${HERWIG_API}")


#--- searching for Boost (for ThePEG filesystem handling)
find_package(Boost 1.62 COMPONENTS system filesystem)
if(NOT Boost_FOUND)
  message(STATUS "Herwig was found but the version of Boost is incompatible.")
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources HerwigHadroniser.cpp HerwigEventInterface.cpp)

#----- build the object

cepgen_build(CepGenHerwig SOURCES ${sources}
    EXT_LIBS ${Boost_LIBRARIES} ${THEPEG} ${THEPEG_LHE} ${HERWIG} ${HERWIG_API}
    EXT_HEADERS ${THEPEG_INCLUDE} ${HERWIG_INCLUDE} ${HEPMC_INCLUDE}
    PROPERTY CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation"
    INSTALL_COMPONENT herwig)

#----- extra installation rules

install(DIRECTORY DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT herwig FILES_MATCHING PATTERN "*.h")
