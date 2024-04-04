set(EPIC_DIRS $ENV{EPIC_DIR})
find_path(EPIC_SOURCES_DIR NAMES Epic.cpp HINTS ${EPIC_DIRS} PATH_SUFFIXES src)
find_path(EPIC_INCLUDE NAMES Epic.h HINTS ${EPIC_DIRS} PATH_SUFFIXES include)

if(NOT EPIC_SOURCES_DIR OR NOT EPIC_INCLUDE)
  return()
endif()

find_path(PARTONS_DIR NAMES include/partons/Partons.h HINTS $ENV{PARTONS_DIR} PATH_SUFFIXES partons)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${PARTONS_DIR}/cmake/Modules ${PARTONS_DIR}/partons/cmake/Modules)
find_package(PARTONS)
if(NOT PARTONS_FOUND)
  message(STATUS "EpIC was found, but not PARTONS.")
  return()
endif()
find_package(GSL)
find_package(SFML COMPONENTS system REQUIRED)
find_package(ElementaryUtils REQUIRED)
find_package(NumA++ REQUIRED)
find_package(Apfel++ REQUIRED)
find_package(ROOT COMPONENTS Core Hist Minuit Foam REQUIRED)
find_package(HepMC3 REQUIRED)

file(GLOB_RECURSE epic_sources ${EPIC_SOURCES_DIR}/*.cpp)

#----- build the object
cepgen_build(CepGenEpIC SOURCES ${epic_sources}
    EXT_HEADERS ${EPIC_INCLUDE} ${GSL_INCLUDE_DIRS} ${QT_INCLUDE_DIRS} ${SFML_INCLUDE_DIR} ${ElementaryUtils_INCLUDE_DIR} ${NumA++_INCLUDE_DIR} ${PARTONS_INCLUDE_DIR} ${Apfel++_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${HEPMC3_INCLUDE_DIR}
    EXT_LIBS ${GSL_LIBRARIES} ${SFML_LIBRARIES} ${ElementaryUtils_LIBRARIES} ${NumA++_LIBRARIES} ${PARTONS_LIBRARIES} ${Apfel++_LIBRARIES} ${ROOT_LIBRARIES} ${HEPMC3_LIBRARIES}
    OPTIONS "-Wno-deprecated-copy"
    INSTALL_COMPONENT epic)
