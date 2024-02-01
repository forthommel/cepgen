set(PARTONS_DIRS $ENV{PARTONS_DIR} $ENV{PARTONS} /usr /usr/local)
find_library(PARTONS_LIB NAMES PARTONS HINTS ${PARTONS_DIRS} PATH_SUFFIXES lib)
find_path(PARTONS_INCLUDE NAMES partons/Partons.h HINTS ${PARTONS_DIRS} PATH_SUFFIXES include)

if(NOT PARTONS_LIB)
  return()
endif()

#----- build the object

cepgen_build(CepGenPartons SOURCES *.cpp
    EXT_LIBS ${PARTONS_LIB}
    EXT_HEADERS ${PARTONS_INCLUDE}
    INSTALL_COMPONENT partons)
cpack_add_component(partons
    DISPLAY_NAME "CepGen Partons wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Partons library"
    DEPENDS lib)
