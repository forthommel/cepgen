#--- searching for MARTY
find_library(MARTY marty HINTS $ENV{MARTY_DIR} PATH_SUFFIXES lib build)
find_path(MARTY_INCLUDE marty.h HINTS $ENV{MARTY_DIR} PATH_SUFFIXES src)

if(NOT MARTY)
  return()
endif()

#----- build the object

cepgen_build(CepGenMarty SOURCES *.cpp
    EXT_LIBS ${MARTY}
    EXT_HEADERS ${MARTY_INCLUDE}
    PROPERTY CXX_STANDARD 17
    INSTALL_COMPONENT marty)
