#--- searching for HepLibF
find_library(HEPLIB HepLib HINTS $ENV{HEPLIB_PATH} PATH_SUFFIXES lib)
find_path(HEPLIB_INCLUDE HepLib.h HINTS $ENV{HEPLIB_PATH} PATH_SUFFIXES include)

if(NOT HEPLIB)
  return()
endif()

#----- build the object

cepgen_build(CepGenHepLib
    SOURCES *.cpp
    TESTS test/*.cc
    EXT_LIBS ${HEPLIB}
    EXT_HEADERS ${HEPLIB_INCLUDE}
    PROPERTY CXX_STANDARD 17
    CFLAGS "-Wno-pedantic -Wno-unused-parameter -Wno-sign-compare"
    INSTALL_COMPONENT heplib)
