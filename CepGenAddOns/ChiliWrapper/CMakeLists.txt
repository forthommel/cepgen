set(CHILI_DIRS $ENV{CHILI_DIR} /usr /usr/local)
find_library(CHILI NAMES Chili HINTS ${CHILI_DIRS} PATH_SUFFIXES lib lib64 build/lib)
find_path(CHILI_INCLUDE NAMES Chili/Channel/Channel.hh HINTS ${CHILI_DIRS} PATH_SUFFIXES include)

if(NOT CHILI AND NOT CHILI_INCLUDE)
  return()
endif()

#----- build the object

cepgen_build(CepGenChili SOURCES *.cpp
    EXT_LIBS ${CHILI}
    EXT_HEADERS ${CHILI_INCLUDE}
    INSTALL_COMPONENT chili)
