find_package(pybind11)
if(NOT pybind11_FOUND)
  return()
endif()

#----- build the object

cepgen_build(CepGenPyBind11 SOURCES *.cpp
    EXT_LIBS ${pybind11_LIBRARIES}
    EXT_HEADERS ${pybind11_INCLUDE_DIRS}
    INSTALL_COMPONENT python)
