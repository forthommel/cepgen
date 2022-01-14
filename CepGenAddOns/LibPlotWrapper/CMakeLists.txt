#--- linking with libplot/plotutils

set(LIBPLOT_DIRS $ENV{LIBPLOT_DIR} /usr /usr/local)
find_library(LIBPLOT NAMES libplotter HINTS ${LIBPLOT_DIRS} PATH_SUFFIXES lib lib64)
find_path(LIBPLOT_INCLUDE NAMES plotter.h HINTS ${LIBPLOT_DIRS} PATH_SUFFIXES src)

message(STATUS ":::${LIBPLOT}")
if(NOT LIBPLOT)
  return()
endif()

#----- build the object

cepgen_build(CepGenLibPlot SOURCES *.cpp
    EXT_LIBS ${LIBPLOT}
    EXT_HEADERS ${LIBPLOT_INCLUDE}
    INSTALL_COMPONENT libplot)
