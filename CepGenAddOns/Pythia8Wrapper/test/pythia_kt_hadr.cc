#include <Pythia8/Pythia.h>

#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Common/EventUtils.h"
#include "CepGenAddOns/Pythia8Wrapper/CepGenEventInterface.h"

using namespace std;

int main(int argc, char* argv[]) {
  auto argparse = cepgen::ArgumentsParser(argc, argv).parse();

  auto gen = cepgen::Generator();

  auto evt = cepgen::utils::generateLPAIREvent();

  auto py_evt = Pythia8::CepGenEventInterface();
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

  CG_DEBUG("main") << evt;

  for (const auto& mode :
       {Pythia8::CepGenEventInterface::Type::central | Pythia8::CepGenEventInterface::Type::partonsKT,
        Pythia8::CepGenEventInterface::Type::central | Pythia8::CepGenEventInterface::Type::partonsKT |
            Pythia8::CepGenEventInterface::Type::beamRemnants,
        Pythia8::CepGenEventInterface::Type::central | Pythia8::CepGenEventInterface::Type::partonsCollinear}) {
    py_evt.feedEvent(evt, mode);
    py_evt.listEvent();
    CG_DEBUG("main").log([&py_evt](auto& log) { py_evt.dumpCorresp(log.stream()); });

    for (int i = 0; i < py_evt.sizePart(); ++i) {
      if (py_evt.cepgenId(i) == Pythia8::CepGenEventInterface::INVALID_ID)
        continue;
      const auto moth1 = *evt[py_evt.cepgenId(i)].mothers().begin(),
                 moth2 = *evt[py_evt.cepgenId(i)].mothers().rbegin();
      CG_TEST_EQUAL(
          py_evt.mother1(i), py_evt.pythiaId(moth1), cepgen::utils::format("[mode %d] part.%d mother 1", mode, i));
      CG_TEST_EQUAL(
          py_evt.mother2(i), py_evt.pythiaId(moth2), cepgen::utils::format("[mode %d] part.%d mother 2", mode, i));
    }
  }

  gen.parametersRef().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  cepgen::ParametersList pythia8_params;
  if (argparse.debugging())
    pythia8_params.set<bool>("debug", true);

  auto cg_pythia = cepgen::EventModifierFactory::get().build("pythia8", pythia8_params);
  cg_pythia->readStrings({"BeamRemnants:primordialKT = off"});
  cg_pythia->setCrossSection(cepgen::Value{1.46161e-1, 1.25691e-3});
  cg_pythia->initialise(gen.parametersRef());
  double evt_weight = 1.;
  cg_pythia->run(evt, evt_weight, true);

  CG_DEBUG("main") << "Pythia 8-filtered event:\n" << evt;

  CG_TEST_EQUAL(evt_weight, 1., "event weight");
  CG_TEST(evt(cepgen::Particle::Role::OutgoingBeam1).size() > 1, "decayed diffractive beam system");
  CG_TEST(evt(cepgen::Particle::Role::OutgoingBeam2).size() == 1, "undecayed elastic beam system");
  cepgen::Momentum daugh_total_momentum;
  for (const auto& daugh : evt.stableDaughters(evt(cepgen::Particle::Role::OutgoingBeam1)[0], true))
    daugh_total_momentum += daugh.get().momentum();
  CG_TEST_EQUIV((daugh_total_momentum - evt(cepgen::Particle::Role::OutgoingBeam1)[0].momentum()).p(),
                0.,
                "diffractive system momentum balance");

  CG_TEST_SUMMARY;
}
