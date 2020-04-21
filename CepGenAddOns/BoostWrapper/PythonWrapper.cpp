#include "CepGen/Version.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
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
      else
        cepgen::proc::Process::prepareKinematics();
    }
    void fillKinematics( bool sym ) override {
      if ( py::override ov = this->get_override( "fillKinematics" ) )
        ov( sym );
      else
        cepgen::proc::Process::fillKinematics( sym );
    }
  };

  struct EventModifierWrap : cepgen::EventModifier, py::wrapper<cepgen::EventModifier>
  {
    void init() override {
      if ( py::override ov = this->get_override( "init" ) )
        ov();
      else
        cepgen::EventModifier::init();
    }
    bool run( cepgen::Event& ev, double& weight, bool full ) override {
      if ( py::override ov = this->get_override( "run" ) )
        return ov( ev, weight, full );
      return cepgen::EventModifier::run( ev, weight, full );
    }
  };

  struct ExportModuleWrap : cepgen::io::ExportModule, py::wrapper<cepgen::io::ExportModule>
  {
    void initialise( const cepgen::Parameters& params ) override {
      if ( py::override ov = this->get_override( "initialise" ) )
        ov( params );
      else
        cepgen::io::ExportModule::initialise( params );
    }
    void operator<<( const cepgen::Event& evt ) override {
      if ( py::override ov = this->get_override( "operator<<" ) )
        ov( evt );
      else
        cepgen::io::ExportModule::operator<<( evt );
    }
  };
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( pl_get_int_ov, cepgen::ParametersList::get<int>, 1, 2 )

BOOST_PYTHON_MODULE( pycepgen )
{
  //----- Process implementation

  py::class_<ProcessWrap, std::shared_ptr<cepgen::proc::Process>, boost::noncopyable>( "Process", "Physics process implementation", py::no_init )
    .def( "computeWeight", py::pure_virtual( &cepgen::proc::Process::computeWeight ) )
    .def( "prepareKinematics", py::pure_virtual( &cepgen::proc::Process::prepareKinematics ) )
    .def( "fillKinematics", py::pure_virtual( &cepgen::proc::Process::fillKinematics ) )
  ;

  //----- Event content

  py::enum_<cepgen::Particle::Status>( "Status", "Particle current status" )
    .value( "primordialIncoming", cepgen::Particle::Status::PrimordialIncoming )
    .value( "debugResonance", cepgen::Particle::Status::DebugResonance )
    .value( "resonance", cepgen::Particle::Status::Resonance )
    .value( "fragmented", cepgen::Particle::Status::Fragmented )
    .value( "propagator", cepgen::Particle::Status::Propagator )
    .value( "incoming", cepgen::Particle::Status::Incoming )
    .value( "undefined", cepgen::Particle::Status::Undefined )
    .value( "finalState", cepgen::Particle::Status::FinalState )
    .value( "undecayed", cepgen::Particle::Status::Undecayed )
    .value( "unfragmented", cepgen::Particle::Status::Unfragmented )
  ;

  py::enum_<cepgen::Particle::Role>( "Role", "Particle role in event" )
    .value( "unknownRole", cepgen::Particle::Role::UnknownRole )
    .value( "incomingBeam1", cepgen::Particle::Role::IncomingBeam1 )
    .value( "incomingBeam2", cepgen::Particle::Role::IncomingBeam2 )
    .value( "outgoingBeam1", cepgen::Particle::Role::OutgoingBeam1 )
    .value( "outgoingBeam2", cepgen::Particle::Role::OutgoingBeam2 )
    .value( "centralSystem", cepgen::Particle::Role::CentralSystem )
    .value( "intermediate", cepgen::Particle::Role::Intermediate )
    .value( "parton1", cepgen::Particle::Role::Parton1 )
    .value( "parton2", cepgen::Particle::Role::Parton2 )
  ;

  cepgen::Particle& ( cepgen::Particle::*particle_set_status )( cepgen::Particle::Status ) = &cepgen::Particle::setStatus;
  cepgen::Particle& ( cepgen::Particle::*particle_set_pdgid )( cepgen::pdgid_t, short ) = &cepgen::Particle::setPdgId;

  py::class_<cepgen::Particle>( "Particle", "Particle kinematics and general information" )
    .add_property( "id", &cepgen::Particle::id, py::make_function( &cepgen::Particle::setId, py::return_self<>() ), "Particle identifier in event" )
    .add_property( "status", &cepgen::Particle::status, py::make_function( particle_set_status, py::return_self<>() ), "Particle status (stable, decayed, ...)" )
    .add_property( "role", &cepgen::Particle::role, py::make_function( &cepgen::Particle::setRole, py::return_self<>() ), "Particle role in event" )
    .add_property( "pdgId", &cepgen::Particle::pdgId, py::make_function( particle_set_pdgid, py::return_self<>() ), "PDG identifier" )
  ;

  py::class_<cepgen::Event>( "Event", "Event content" )
    .def( py::init<bool>() )
  ;

  //----- Event modifiers

  py::class_<EventModifierWrap, std::shared_ptr<cepgen::EventModifier>, boost::noncopyable>( "EventModifier", "Event modification algorithm", py::no_init )
    .def( "init", py::pure_virtual( &cepgen::EventModifier::init ), "Initialise the algorithm" )
    .def( "run", py::pure_virtual( &cepgen::EventModifier::run ), "Launch the event modification algorithm on one event" )
    .def( "setSeed", &cepgen::EventModifier::setSeed, "Set the random number generator seed if any" )
  ;

  //----- Export modules

  py::class_<ExportModuleWrap, std::shared_ptr<cepgen::io::ExportModule>, boost::noncopyable>( "ExportModule", "Export module", py::no_init )
    .add_property( "name", py::make_function( &cepgen::io::ExportModule::name, py::return_value_policy<py::copy_const_reference>() ), "Module name" )
    .def( "initialise", py::pure_virtual( &cepgen::io::ExportModule::initialise ), "Initialise the module from runtime parameters" )
    .def( "setCrossSection", &cepgen::io::ExportModule::setCrossSection, "Specify the cross section + uncertainty for this run" )
    .def( "setEventNumber", &cepgen::io::ExportModule::setEventNumber, "Specify the event number in this run" )
    .def( "feed", py::pure_virtual( &cepgen::io::ExportModule::operator<< ), "Feed the export module with one event" )
  ;

  //----- Runtime parameters

  cepgen::proc::Process& ( cepgen::Parameters::*process_get )() = &cepgen::Parameters::process;
  void ( cepgen::Parameters::*process_set )( cepgen::proc::Process* ) = &cepgen::Parameters::setProcess;
  cepgen::EventModifiersSequence& ( cepgen::Parameters::*evt_mod_seq )()  = &cepgen::Parameters::eventModifiersSequence;
  void ( cepgen::Parameters::*evt_mod_add )( cepgen::EventModifier* ) = &cepgen::Parameters::addModifier;
  cepgen::ExportModulesSequence& ( cepgen::Parameters::*out_mod_seq )()  = &cepgen::Parameters::outputModulesSequence;
  void ( cepgen::Parameters::*out_mod_add )( cepgen::io::ExportModule* ) = &cepgen::Parameters::addOutputModule;

  py::class_<cepgen::Parameters>( "Parameters", "Collection of runtime parameters for CepGen" )
    //.def( py::self_ns::str( py::self_ns::self ) )
    .add_property( "process", py::make_function( process_get, py::return_internal_reference<>() ), process_set, "Process to generate" )
    .def( "eventModifiers", py::make_function( evt_mod_seq, py::return_internal_reference<>() ), "List of event modification algorithms registered" )
    .def( "addEventModifier", evt_mod_add, "Add a new event modification algorithm reference to the stack" )
    .def( "outputModules", py::make_function( out_mod_seq, py::return_internal_reference<>() ), "List of output modules registered" )
    .def( "addOutputModule", out_mod_add, "Add a new output module reference to the stack" )
  ;

  //----- Utilities

  py::class_<cepgen::ParametersList>( "ParametersList", "Structured handler for mixed parameters" )
    .add_property( "keys", &cepgen::ParametersList::keys, "List of all parameters indexes" )
    //.def( "get", &cepgen::ParametersList::get<int>, pl_get_int_ov(), "Retrieve a parameter by its index" )
    //.def( "get", &( cepgen::ParametersList::get<int> ), pl_get_int_ov(), "Retrieve a parameter by its index" )
    //.def( "get", cepgen::ParametersList::get<int>, pl_get_int_ov(), "Retrieve a parameter by its index" )
    //.def( "get", get<int>, pl_get_int_ov(), "Retrieve a parameter by its index" )
  ;
}
