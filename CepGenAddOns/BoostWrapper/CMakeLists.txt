find_package(Boost COMPONENTS python)

if(NOT Boost_FOUND)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

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

find_package(PythonLibs 2.7 REQUIRED)
PYTHON_ADD_MODULE(pycepgen PythonWrapper.cpp)
target_include_directories(pycepgen PRIVATE ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
target_link_libraries(pycepgen CepGen ${Boost_LIBRARIES})
install(TARGETS pycepgen DESTINATION ${PYTHON_INSTALL_PATH})
#--- copy the __init__.py
#configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py COPYONLY)
#install(TARGETS pycepgen __init__.py DESTINATION ${PYTHON_INSTALL_PATH})

#----- build the object

cepgen_build(CepGenBoost SOURCES ${sources}
    EXT_LIBS ${Boost_LIBRARIES}
    EXT_HEADERS ${Boost_INCLUDE_DIRS}
    CFLAGS "-pthread -Wno-pedantic"
    INSTALL_COMPONENT boost)
