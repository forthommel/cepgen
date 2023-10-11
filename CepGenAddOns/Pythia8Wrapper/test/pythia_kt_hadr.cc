#include <Pythia8/Pythia.h>

#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Common/EventUtils.h"
#include "CepGenAddOns/Pythia8Wrapper/PythiaEventInterface.h"

using namespace std;

int main() {
  //cepgen::initialise();
  auto gen = cepgen::Generator();

  auto evt = cepgen::utils::generateLPAIREvent();

  auto py_evt = Pythia8::CepGenEvent();
  auto pos_beam = cepgen::Beam({})
                      .setPdgId(evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1).pdgId())
                      .setMomentum(evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1).momentum())
                      .setElastic(false);
  auto neg_beam = cepgen::Beam({})
                      .setPdgId(evt.oneWithRole(cepgen::Particle::Role::IncomingBeam2).pdgId())
                      .setMomentum(evt.oneWithRole(cepgen::Particle::Role::IncomingBeam2).momentum())
                      .setElastic(true);
  py_evt.initialise(pos_beam, neg_beam);

  CG_TEST_EQUAL((unsigned long long)py_evt.idBeamA(),
                evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1).pdgId(),
                "first incoming beam PDG id");
  CG_TEST_EQUAL((unsigned long long)py_evt.idBeamB(),
                evt.oneWithRole(cepgen::Particle::Role::IncomingBeam2).pdgId(),
                "second incoming beam PDG id");
  CG_TEST_EQUAL(py_evt.eBeamA(),
                evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1).momentum().pz(),
                "first incoming beam pz");
  CG_TEST_EQUAL(py_evt.eBeamB(),
                evt.oneWithRole(cepgen::Particle::Role::IncomingBeam2).momentum().pz(),
                "second incoming beam pz");

  evt.dump();
  for (const auto& mode : {Pythia8::CepGenEvent::Type::centralAndPartons,
                           Pythia8::CepGenEvent::Type::centralAndFullBeamRemnants,
                           Pythia8::CepGenEvent::Type::centralAndBeamRemnants}) {
    py_evt.feedEvent(evt, mode);
    py_evt.listEvent();
    py_evt.dumpCorresp();
    for (int i = 0; i < py_evt.sizePart(); ++i) {
      if (py_evt.cepgenId(i) == Pythia8::CepGenEvent::INVALID_ID)
        continue;
      const auto moth1 = *evt[py_evt.cepgenId(i)].mothers().begin(),
                 moth2 = *evt[py_evt.cepgenId(i)].mothers().rbegin();
      CG_LOG << moth1 << "(" << py_evt.pythiaId(moth2) << "):" << moth2 << "(" << py_evt.pythiaId(moth1) << ")";
      CG_TEST_EQUAL(
          py_evt.mother1(i), py_evt.pythiaId(moth1), cepgen::utils::format("[mode %d] part.%d mother 1", (int)mode, i));
      CG_TEST_EQUAL(
          py_evt.mother2(i), py_evt.pythiaId(moth2), cepgen::utils::format("[mod3 %d] part.%d mother 2", (int)mode, i));
    }
  }

  gen.parametersRef().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set("kinematics",
                                   cepgen::ParametersList().set<double>("cmEnergy", 13.e3).set<int>("mode", 3))));

  auto cg_pythia = cepgen::EventModifierFactory::get().build("pythia8");
  cg_pythia->setCrossSection(cepgen::Value{1., 1.e-3});
  cg_pythia->initialise(gen.parametersRef());
  double evt_weight = 1.;
  cg_pythia->run(evt, evt_weight, true);

  CG_TEST_SUMMARY;
}
