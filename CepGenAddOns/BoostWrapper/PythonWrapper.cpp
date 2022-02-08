#include <boost/config.hpp>
#include <boost/python.hpp>
#include <boost/version.hpp>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Version.h"

namespace {
  namespace py = boost::python;

  struct ProcessWrap : cepgen::proc::Process, py::wrapper<cepgen::proc::Process> {
    explicit ProcessWrap() : cepgen::proc::Process(cepgen::ParametersList()) {}
    double computeWeight() override {
      if (py::override ov = this->get_override("computeWeight"))
        return ov();
      return cepgen::proc::Process::computeWeight();
    }
    void prepareKinematics() override {
      if (py::override ov = this->get_override("prepareKinematics"))
        ov();
      else
        cepgen::proc::Process::prepareKinematics();
    }
    void fillKinematics(bool sym) override {
      if (py::override ov = this->get_override("fillKinematics"))
        ov(sym);
      else
        cepgen::proc::Process::fillKinematics(sym);
    }
  };

  struct EventModifierWrap : cepgen::EventModifier, py::wrapper<cepgen::EventModifier> {
    void initialise() override {
      if (py::override ov = this->get_override("initialise"))
        ov();
      else
        cepgen::EventModifier::initialise();
    }
    bool run(cepgen::Event& ev, double& weight, bool full) override {
      if (py::override ov = this->get_override("run"))
        return ov(ev, weight, full);
      return cepgen::EventModifier::run(ev, weight, full);
    }
  };

  struct EventExporterWrap : cepgen::EventExporter, py::wrapper<cepgen::EventExporter> {
    void initialise() override {
      if (py::override ov = this->get_override("initialise"))
        ov();
      else
        cepgen::EventExporter::initialise();
    }
    void operator<<(const cepgen::Event& evt) override {
      if (py::override ov = this->get_override("operator<<"))
        ov(evt);
      else
        cepgen::EventExporter::operator<<(evt);
    }
  };

  PyObject* except_type = nullptr;

  void translate_exception(const cepgen::Exception& e) {
    if (except_type == NULL)
      return;
    PyErr_SetObject(except_type, py::object(e).ptr());
  }

  /*std::shared_ptr<cepgen::ModuleFactory> getModuleFactory() {
    return std::shared_ptr<cepgen::ModuleFactory>(&cepgen::ModuleFactory::get(), [](const void*) {});
  }*/

  std::shared_ptr<cepgen::StructureFunctionsFactory> getStructureFunctionsFactory() {
    return std::shared_ptr<cepgen::StructureFunctionsFactory>(&cepgen::StructureFunctionsFactory::get(),
                                                              [](const void*) {});
  }
}  // namespace

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(part_set_mom_ov, cepgen::Particle::setMomentum, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mod_build_ov, build, 1, 2)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(pl_get_int_ov, get<int>, 1, 2)

BOOST_PYTHON_MODULE(pycepgen) {
  //----- Process implementation

  py::class_<ProcessWrap, std::shared_ptr<cepgen::proc::Process>, boost::noncopyable>(
      "Process", "Physics process implementation", py::no_init)
      .def("computeWeight", py::pure_virtual(&cepgen::proc::Process::computeWeight))
      .def("prepareKinematics", py::pure_virtual(&cepgen::proc::Process::prepareKinematics))
      .def("fillKinematics", py::pure_virtual(&cepgen::proc::Process::fillKinematics));

  //----- Event content

  py::enum_<cepgen::Particle::Status>("ParticleStatus", "Particle current status")
      .value("primordialIncoming", cepgen::Particle::Status::PrimordialIncoming)
      .value("debugResonance", cepgen::Particle::Status::DebugResonance)
      .value("resonance", cepgen::Particle::Status::Resonance)
      .value("fragmented", cepgen::Particle::Status::Fragmented)
      .value("propagator", cepgen::Particle::Status::Propagator)
      .value("incoming", cepgen::Particle::Status::Incoming)
      .value("undefined", cepgen::Particle::Status::Undefined)
      .value("finalState", cepgen::Particle::Status::FinalState)
      .value("undecayed", cepgen::Particle::Status::Undecayed)
      .value("unfragmented", cepgen::Particle::Status::Unfragmented);

  py::enum_<cepgen::Particle::Role>("ParticleRole", "Particle role in event")
      .value("unknownRole", cepgen::Particle::Role::UnknownRole)
      .value("incomingBeam1", cepgen::Particle::Role::IncomingBeam1)
      .value("incomingBeam2", cepgen::Particle::Role::IncomingBeam2)
      .value("outgoingBeam1", cepgen::Particle::Role::OutgoingBeam1)
      .value("outgoingBeam2", cepgen::Particle::Role::OutgoingBeam2)
      .value("centralSystem", cepgen::Particle::Role::CentralSystem)
      .value("intermediate", cepgen::Particle::Role::Intermediate)
      .value("parton1", cepgen::Particle::Role::Parton1)
      .value("parton2", cepgen::Particle::Role::Parton2);

  py::class_<cepgen::Particle>("Particle", "Particle kinematics and general information")
      .def(py::self_ns::str(py::self_ns::self))
      .add_property("id",
                    &cepgen::Particle::id,
                    py::make_function(&cepgen::Particle::setId, py::return_self<>()),
                    "Particle identifier in event")
      .add_property("status",
                    &cepgen::Particle::status,
                    py::make_function(static_cast<cepgen::Particle& (cepgen::Particle::*)(cepgen::Particle::Status)>(
                                          &cepgen::Particle::setStatus),
                                      py::return_self<>()),
                    "Particle status (stable, decayed, ...)")
      .add_property("role",
                    &cepgen::Particle::role,
                    py::make_function(&cepgen::Particle::setRole, py::return_self<>()),
                    "Particle role in event")
      .add_property("pdgId",
                    &cepgen::Particle::pdgId,
                    py::make_function(static_cast<cepgen::Particle& (cepgen::Particle::*)(cepgen::pdgid_t, short)>(
                                          &cepgen::Particle::setPdgId),
                                      py::return_self<>()),
                    "PDG identifier")
      .add_property(
          "momentum",
          py::make_function(static_cast<cepgen::Momentum& (cepgen::Particle::*)()>(&cepgen::Particle::momentum),
                            py::return_internal_reference<>()),
          +[](cepgen::Particle& part, const cepgen::Momentum& mom) { part.setMomentum(mom, true); },
          "4-momentum");

  py::class_<cepgen::Event>("Event", "Event content").def(py::init<bool>()).def(py::self_ns::str(py::self_ns::self))
      //.def( py::vector_indexing_suite<cepgen::Particles>() )
      //.def("addParticle",
      //     py::make_function(
      //         static_cast<cepgen::Particle& (cepgen::Event::*)(cepgen::Particle&, bool)>(&cepgen::Event::addParticle),
      //         py::return_value_policy<py::reference_existing_object>()),
      //     "Add a particle object reference to the event")
      ;

  //----- Event modifiers

  py::class_<EventModifierWrap, std::shared_ptr<cepgen::EventModifier>, boost::noncopyable>(
      "EventModifier", "Event modification algorithm", py::no_init)
      //.def("initialise", py::pure_virtual(&cepgen::EventModifier::initialise), "Initialise the algorithm")
      .def("initialise", &cepgen::EventModifier::initialise, "Initialise the algorithm")
      .def("run", py::pure_virtual(&cepgen::EventModifier::run), "Launch the event modification algorithm on one event")
      .def("setSeed", &cepgen::EventModifier::setSeed, "Set the random number generator seed if any");

  //----- Export modules

  py::class_<EventExporterWrap, std::shared_ptr<cepgen::EventExporter>, boost::noncopyable>(
      "EventExporter", "Export module", py::no_init)
      .add_property(
          "name",
          py::make_function(&cepgen::EventExporter::name, py::return_value_policy<py::copy_const_reference>()),
          "Module name")
      .def("initialise",
           //py::pure_virtual(&cepgen::EventExporter::initialise),
           &cepgen::EventExporter::initialise)  //,
                                                //"Initialise the module from runtime parameters")
      .def("setCrossSection",
           &cepgen::EventExporter::setCrossSection,
           "Specify the cross section + uncertainty for this run")
      .def("setEventNumber", &cepgen::EventExporter::setEventNumber, "Specify the event number in this run")
      .def("feed", py::pure_virtual(&cepgen::EventExporter::operator<<), "Feed the export module with one event");

  //----- Runtime parameters

  py::class_<cepgen::Parameters>("Parameters", "Collection of runtime parameters for CepGen")
      //.def(py::self_ns::str(py::self_ns::self))
      .def("eventModifiers",
           py::make_function(static_cast<cepgen::EventModifiersSequence& (cepgen::Parameters::*)()>(
                                 &cepgen::Parameters::eventModifiersSequence),
                             py::return_internal_reference<>()),
           "List of event modification algorithms registered")
      .def("addEventModifier",
           static_cast<void (cepgen::Parameters::*)(cepgen::EventModifier*)>(&cepgen::Parameters::addModifier),
           "Add a new event modification algorithm reference to the stack")
      .def("eventExporters",
           py::make_function(static_cast<cepgen::EventExportersSequence& (cepgen::Parameters::*)()>(
                                 &cepgen::Parameters::eventExportersSequence),
                             py::return_internal_reference<>()),
           "List of output modules registered")
      .def("addEventExporter",
           static_cast<void (cepgen::Parameters::*)(cepgen::EventExporter*)>(&cepgen::Parameters::addEventExporter),
           "Add a new output module reference to the stack")
      .add_property(
          "process",
          py::make_function(static_cast<cepgen::proc::Process& (cepgen::Parameters::*)()>(&cepgen::Parameters::process),
                            py::return_internal_reference<>()),
          static_cast<void (cepgen::Parameters::*)(cepgen::proc::Process*)>(&cepgen::Parameters::setProcess),
          "Process to generate");

  //----- Structure functions

  py::class_<cepgen::strfun::Parameterisation>("StructureFunctions", "Structure functions parameterisation")
      .def("F1", &cepgen::strfun::Parameterisation::F1)
      .def("F2", &cepgen::strfun::Parameterisation::F2)
      .def("FL", &cepgen::strfun::Parameterisation::FL)
      .def("W1", &cepgen::strfun::Parameterisation::W1)
      .def("W2", &cepgen::strfun::Parameterisation::W2)
      .def("FE", &cepgen::strfun::Parameterisation::FE)
      .def("FM", &cepgen::strfun::Parameterisation::FM);
  py::register_ptr_to_python<cepgen::strfun::Parameterisation*>();

  py::class_<cepgen::StructureFunctionsFactory, boost::noncopyable>("StructureFunctionsFactory", py::no_init)
      .def(
          "build", +[](int mod) { return cepgen::StructureFunctionsFactory::get().build(mod); })
      .def(
          "build",
          +[](const cepgen::ParametersList& plist) { return cepgen::StructureFunctionsFactory::get().build(plist); });
  //static_cast<std::unique_ptr<cepgen::strfun::Parameterisation> (cepgen::StructureFunctionsFactory::*)(
  //         const int&, const cepgen::ParametersList&)>(&cepgen::StructureFunctionsFactory::build));

  //----- Utilities

  py::class_<cepgen::Momentum>("Momentum", "4-momentum container")
      .def(py::init<double, double, double, double>())
      .def(py::self * int())
      .def(py::self * double())
      .def(py::self + py::self)
      .def(py::self - py::self)
      //.def(py::self == py::self)
      .def(py::self_ns::str(py::self_ns::self))
      .add_property("px",
                    &cepgen::Momentum::px,
                    py::make_function(&cepgen::Momentum::setPx, py::return_self<>()),
                    "Horizontal momentum (GeV/c)")
      .add_property("py",
                    &cepgen::Momentum::py,
                    py::make_function(&cepgen::Momentum::setPy, py::return_self<>()),
                    "Vertical momentum (GeV/c)")
      .add_property("pz",
                    &cepgen::Momentum::pz,
                    py::make_function(&cepgen::Momentum::setPz, py::return_self<>()),
                    "Longitudinal momentum (GeV/c)")
      .add_property("p", &cepgen::Momentum::p, "3-momentum norm (GeV/c)")
      .add_property("pt", &cepgen::Momentum::pt, "Transverse momentum (GeV/c)")
      .add_property("eta", &cepgen::Momentum::eta, "Pseudorapidity")
      .add_property("rapidity", &cepgen::Momentum::rapidity, "Rapidity")
      .add_property("phi", &cepgen::Momentum::phi, "Azimuthal angle (rad)")
      .add_property("energy",
                    &cepgen::Momentum::energy,
                    py::make_function(&cepgen::Momentum::setEnergy, py::return_self<>()),
                    "Energy (GeV)")
      .add_property("mass", &cepgen::Momentum::mass, "Mass (GeV/c^2)");

  py::class_<cepgen::ParametersList>("ParametersList", "Structured handler for mixed parameters")
      .def(py::self_ns::str(py::self_ns::self))
      .add_property("keys", &cepgen::ParametersList::keys, "List of all parameters indexes")
      /*.def("get", &cepgen::ParametersList::get<int>, pl_get_int_ov(), "Retrieve a parameter by its index")
      .def("get", &(cepgen::ParametersList::get<int>), pl_get_int_ov(), "Retrieve a parameter by its index")
      .def("get", cepgen::ParametersList::get<int>, pl_get_int_ov(), "Retrieve a parameter by its index")
      .def("get", get<int>, pl_get_int_ov(), "Retrieve a parameter by its index")*/
      ;

  py::class_<cepgen::Exception, boost::noncopyable> except("Exception", "Generic exception", py::no_init);
  except.add_property("message", &cepgen::Exception::message, "Human-readable error message");
  except_type = except.ptr();
  py::register_exception_translator<cepgen::Exception>(&translate_exception);
}
