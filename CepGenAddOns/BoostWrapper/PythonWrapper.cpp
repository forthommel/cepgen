#include "CepGen/Version.h"
#include "CepGen/Parameters.h"

#include <boost/version.hpp>
#include <boost/config.hpp>
#include <boost/python.hpp>

namespace
{
  namespace py = boost::python;
}

BOOST_PYTHON_MODULE( pycepgen )
{
  py::class_<cepgen::Parameters>( "Parameters", "Collection of runtime parameters for CepGen" )
    .def( py::init() )
    .def( py::self_ns::str( py::self_ns::self )
    .add_property( "process", &cepgen::Parameters::process, &cepgen::Parameters::setProcess, "Process to generate" )
  ;
}
