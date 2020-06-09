find_package(Boost)

if(NOT Boost_FOUND)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

<<<<<<< HEAD
set(sources BoostTreeHandler.cpp
            BoostTreeUtils.cpp
            BoostTrapAnalyticalIntegrator.cpp
            BoostGaussKronrodAnalyticalIntegrator.cpp
            BoostGaussLegendreAnalyticalIntegrator.cpp)
file(GLOB boost_int_include ${Boost_INCLUDE_DIRS}/boost/math/quadrature/naive_monte_carlo.hpp)
file(GLOB boost_py_include ${Boost_INCLUDE_DIRS}/boost/python.hpp)
if(boost_int_include)
  list(APPEND sources IntegratorNaive.cpp)
endif()
if(boost_py_include)
  list(APPEND sources PythonConfigWriter.cpp)
endif()

#--- Python API
find_package(PythonLibs)
if(PYTHONLIBS_FOUND)
  find_package(PythonInterp)
  set(BOOST_COMPONENTS system)
  if(PYTHONINTERP_FOUND)
    list(APPEND BOOST_COMPONENTS
         "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
  endif()
  find_package(Boost COMPONENTS system OPTIONAL_COMPONENTS ${BOOST_COMPONENTS})
  PYTHON_ADD_MODULE(pycepgen PythonWrapper.cpp)
  target_include_directories(pycepgen PRIVATE ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(pycepgen CepGen ${Boost_LIBRARIES})
  #configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py COPYONLY)
  execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(TARGETS pycepgen DESTINATION ${PYTHON_INSTALL_PATH} COMPONENT boost)
endif()

#----- build the object

cepgen_build(CepGenBoost SOURCES ${sources}
    EXT_LIBS ${Boost_LIBRARIES}
    EXT_HEADERS ${Boost_INCLUDE_DIRS}
    CFLAGS "-pthread -Wno-pedantic"
    INSTALL_COMPONENT boost)
