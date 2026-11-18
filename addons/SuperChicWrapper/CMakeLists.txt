#--- searching for SuperChic
set(SUPERCHIC_DIRS $ENV{SUPERCHIC_DIR} /usr /usr/local /opt/pythia8)
find_library(SUPERCHIC superchicLib HINTS ${SUPERCHIC_DIRS} PATH_SUFFIXES lib lib64)

if(NOT SUPERCHIC)
  return()
endif()

#----- build the object
cepgen_build(CepGenSuperChic
    SOURCES src/*.cpp
    LIBRARIES ${SUPERCHIC}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    COMPONENT superchic)
cpack_add_component(superchic
    DISPLAY_NAME "CepGen SuperChic wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to SuperChic"
    DEPENDS lib)
#install(DIRECTORY CepGenSuperChic
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#    COMPONENT superchic)
