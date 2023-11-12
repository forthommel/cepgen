#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  auto gen = cepgen::Generator();
  gen.parametersRef().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  auto pythia = cepgen::EventModifierFactory::get().build("pythia8");
  pythia->setCrossSection(cepgen::Value{1., 0.});
  pythia->initialise(gen.parametersRef());

  auto evt = cepgen::Event::minimal(1);
  auto& tau = evt[cepgen::Particle::CentralSystem][0].get();
  tau.setStatus(cepgen::Particle::Status::Undecayed);
  tau.setPdgId(cepgen::PDG::tau);
  tau.setMomentum(cepgen::Momentum(0., 0., 1000.), false);
  for (const auto& role : {cepgen::Particle::Parton1, cepgen::Particle::Parton2}) {
    auto& part = evt[role][0].get();
    part.setPdgId(cepgen::PDG::photon);
    part.setMomentum(tau.momentum() * 0.5);
  }
  const auto evt_size_bef = evt.size();

  double weight;
  pythia->run(evt, weight, true);

  CG_TEST_EQUAL(
      evt(cepgen::Particle::CentralSystem)[0].status(), cepgen::Particle::Status::Resonance, "tau 'decayed' status");
  CG_TEST(evt_size_bef != evt.size(), "decay");

  CG_TEST_SUMMARY;
}
