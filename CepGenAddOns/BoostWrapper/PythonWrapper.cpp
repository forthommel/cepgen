#include "CepGen/Version.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Event/Event.h"

#include <boost/version.hpp>
#include <boost/config.hpp>
#include <boost/python.hpp>

namespace
{
  namespace py = boost::python;

  struct ProcessWrap : cepgen::proc::Process, py::wrapper<cepgen::proc::Process>
  {
    explicit ProcessWrap() :
      cepgen::proc::Process( cepgen::ParametersList() ) {}
    double computeWeight() override {
      if ( py::override ov = this->get_override( "computeWeight" ) )
        return ov();
      return cepgen::proc::Process::computeWeight();
    }
    void prepareKinematics() override {
      if ( py::override ov = this->get_override( "prepareKinematics" ) )
        ov();
      cepgen::proc::Process::prepareKinematics();
    }
    void fillKinematics( bool sym ) override {
      if ( py::override ov = this->get_override( "fillKinematics" ) )
        ov( sym );
      cepgen::proc::Process::fillKinematics( sym );
    }
  };

  struct EventModifierWrap : cepgen::EventModifier, py::wrapper<cepgen::EventModifier>
  {
    void init() override {
      if ( py::override ov = this->get_override( "init" ) )
        ov();
      cepgen::EventModifier::init();
    }
    bool run( cepgen::Event& ev, double& weight, bool full ) override {
      if ( py::override ov = this->get_override( "run" ) )
        return ov( ev, weight, full );
      return cepgen::EventModifier::run( ev, weight, full );
    }
  };
}

BOOST_PYTHON_MODULE( pycepgen )
{
  //----- Process implementation

  py::class_<ProcessWrap, std::shared_ptr<cepgen::proc::Process>, boost::noncopyable>( "Process", "Physics process implementation", py::no_init )
    .def( "computeWeight", py::pure_virtual( &cepgen::proc::Process::computeWeight ) )
    .def( "prepareKinematics", py::pure_virtual( &cepgen::proc::Process::prepareKinematics ) )
    .def( "fillKinematics", py::pure_virtual( &cepgen::proc::Process::fillKinematics ) )
  ;

  //----- Event modifiers

  py::class_<EventModifierWrap, std::shared_ptr<cepgen::EventModifier>, boost::noncopyable>( "EventModifier", "Event modification algorithm", py::no_init )
    .def( "init", py::pure_virtual( &cepgen::EventModifier::init ) )
    .def( "run", py::pure_virtual( &cepgen::EventModifier::run ) )
  ;

  //----- Runtime parameters

  cepgen::proc::Process& ( cepgen::Parameters::*process_get )() = &cepgen::Parameters::process;
  void ( cepgen::Parameters::*process_set )( cepgen::proc::Process* ) = &cepgen::Parameters::setProcess;

  py::class_<cepgen::Parameters>( "Parameters", "Collection of runtime parameters for CepGen" )
    //.def( py::self_ns::str( py::self_ns::self ) )
    .add_property( "process", py::make_function( process_get, py::return_internal_reference<>() ), process_set, "Process to generate" )
  ;
}
