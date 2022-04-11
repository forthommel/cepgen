#--- searching for Herwig 6
set(HERWIG6_DIRS $ENV{HERWIG6_DIR} ${HERWIG6_DIR})
find_library(HERWIG6 herwig HINTS ${HERWIG6_DIRS} PATH_SUFFIXES lib)

set(HERWIG6LIB)
if(HERWIG6)
  list(APPEND HERWIG6LIB ${HERWIG6})
else()
  file(GLOB HERWIG6_SRC ${PROJECT_SOURCE_DIR}/External/herwig6*.f)
  if(EXISTS $ENV{HERWIG6_SRC})
    file(GLOB herwig6_src $ENV{HERWIG6_SRC})
  elseif(EXISTS ${HERWIG6_SRC})
    file(GLOB herwig6_src ${HERWIG6_SRC})
  else()
    return()
  endif()
  message(STATUS "Herwig 6 source found in ${herwig6_src}")
  #set(CMAKE_Fortran_FLAGS_RELEASE "-g -O3")
  add_library(herwig6 SHARED ${herwig6_src})
  install(TARGETS herwig6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT hewrig6)
  set_property(TARGET herwig6 PROPERTY Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  #target_compile_definitions(herwig6 PUBLIC -Wno-maybe-uninitialized -Wno-conversion -g)
  target_compile_options(herwig6 PRIVATE $<$<CONFIG:RELEASE>:-g>)
  list(APPEND HERWIG6LIB herwig6)
endif()

#----- build the object

cepgen_build(CepGenHerwig6
    SOURCES *.cpp
    EXT_LIBS ${HERWIG6LIB}
    INSTALL_COMPONENT herwig6)

