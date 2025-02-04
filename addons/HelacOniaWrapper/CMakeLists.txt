#--- searching for HELAC-Onia
set(HELACOnia_DIRS $ENV{HELACOnia_PATH} $ENV{HELACOnia_DIR} ${HELACOnia_DIR})
find_path(HelacOnia NAMES Helac_Global.f90 HINTS ${HELACOnia_DIRS} PATH_SUFFIXES src)
find_path(HelacOnia_decay NAMES Decay_interface.f90 HINTS ${HELACOnia_DIRS} PATH_SUFFIXES decay)
find_path(HelacOnia_interface NAMES QEDPS_interface.f90 HINTS ${HELACOnia_DIRS} PATH_SUFFIXES shower/interface)
find_path(HelacOnia_plot NAMES plot_user.f90 HINTS ${HELACOnia_DIRS} PATH_SUFFIXES analysis/user)
if (NOT HelacOnia)
  return()
endif ()

file(GLOB HelacOnia_SOURCES ${HelacOnia}/[a-zA-Z0-9]*.f90)
file(GLOB HelacOnia_setscale_default_source ${HelacOnia}/setscale_default.f90)
#file(GLOB HelacOnia_fo_plot_source ${HelacOnia}/FO_plot.f90)
file(GLOB HelacOnia_decay_SOURCES ${HelacOnia_decay}/[a-zA-Z0-9]*.f90)
file(GLOB HelacOnia_qedps_SOURCES ${HelacOnia_interface}/QEDPS_interface.f90)
file(GLOB HelacOnia_plot_SOURCES ${HelacOnia_plot}/plot_user.f90)
#file(GLOB HelacOnia_invalid_sources ${HelacOnia}/._*.f90)
list(REMOVE_ITEM HelacOnia_SOURCES ${HelacOnia_setscale_default_source})
message(STATUS "HELAC-Onia sources: ${HelacOnia_SOURCES}")

#----- build the objects
cepgen_build(CepGenHelacOnia
        SOURCES src/*.cpp ${HelacOnia_SOURCES} ${HelacOnia_decay_SOURCES} ${HelacOnia_qedps_SOURCES} ${HelacOnia_plot_SOURCES}
        OPTIONS "-Wa,-Wno-tabs -Wno-unused-variable -Wno-unused-function -Wno-maybe-uninitialized"
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
        COMPONENT HelacOnia)
cpack_add_component(HelacOnia
        DISPLAY_NAME " CepGen HELAC-Onia wrappers library"
        DESCRIPTION "Collection of CepGen wrappers to the HELAC-Onia library"
        DEPENDS lib)
