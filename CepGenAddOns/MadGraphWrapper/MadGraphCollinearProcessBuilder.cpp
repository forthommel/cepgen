/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcessBuilder.h"

using namespace cepgen;

class MadGraphCollinearProcessBuilder : public MadGraphProcessBuilder, public proc::Process {
public:
  explicit MadGraphCollinearProcessBuilder(const ParametersList&);
  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphCollinearProcessBuilder(*this)); }

  static ParametersDescription description();

  void prepareKinematics() override;
  double computeWeight() override;
  void fillKinematics(bool) override;
  //void prepareProcessKinematics() override;
  //double computeCentralMatrixElement() const override;
  void addEventContent() override;

private:
  std::shared_ptr<collflux::Parameterisation> flux_;
  double x1_{0.}, x2_{0.};
  std::vector<Momentum> momenta_;
};

MadGraphCollinearProcessBuilder::MadGraphCollinearProcessBuilder(const ParametersList& params)
    : MadGraphProcessBuilder(params),
      proc::Process(params, true),
      flux_(collflux::CollinearFluxFactory::get().build(SteeredObject::steer<ParametersList>("collinearFluxes"))) {
  const auto& interm_part = mg5_proc_->intermediatePartons();
  const auto& cent_sys = mg5_proc_->centralSystem();
  CG_DEBUG("MadGraphCollinearProcessBuilder") << "MadGraph_aMC process created for:\n\t"
                                              << "* interm. parts.: " << interm_part << "\n\t"
                                              << "* central system: " << cent_sys << ".";
}

void MadGraphCollinearProcessBuilder::prepareKinematics() {
  defineVariable(x1_, Mapping::linear, {1.e-4, 1.}, {1.e-4, 1.}, "x1");
  defineVariable(x2_, Mapping::linear, {1.e-4, 1.}, {1.e-4, 1.}, "x2");
  if (!mg5_proc_)
    CG_FATAL("MadGraphCollinearProcessBuilder") << "Process not properly linked!";

  CG_INFO("MadGraphCollinearProcessBuilder") << "Preparing process kinematics for card at \"" << params_card_ << "\".";
  mg5_proc_->initialise(params_card_);
}

void MadGraphCollinearProcessBuilder::addEventContent() {
  std::vector<pdgid_t> central_system;
  for (const auto& cs : mg5_proc_->centralSystem())
    central_system.emplace_back(cs);
  proc::Process::setEventContent({{Particle::IncomingBeam1, PDG::proton},
                                  {Particle::IncomingBeam2, PDG::proton},
                                  {Particle::Parton1, mg5_proc_->intermediatePartons()[0]},
                                  {Particle::Parton2, mg5_proc_->intermediatePartons()[1]}},
                                 {{Particle::OutgoingBeam1, {PDG::proton}},
                                  {Particle::OutgoingBeam2, {PDG::proton}},
                                  {Particle::CentralSystem, central_system}});
}

void MadGraphCollinearProcessBuilder::fillKinematics(bool) {
  (*event_)[Particle::Parton1][0].get().setMomentum(momenta_.at(0));
  (*event_)[Particle::Parton2][0].get().setMomentum(momenta_.at(1));
  Momentum central_system;
  for (size_t i = 2; i < momenta_.size(); ++i) {
    (*event_)[Particle::CentralSystem][i - 2].get().setMomentum(momenta_.at(i));
    central_system += momenta_.at(i);
  }
  (*event_)[Particle::Intermediate][0].get().setMomentum(central_system);
  (*event_)[Particle::OutgoingBeam1][0].get().setMomentum((*event_)[Particle::IncomingBeam1][0].get().momentum() -
                                                          momenta_.at(0));
  (*event_)[Particle::OutgoingBeam2][0].get().setMomentum((*event_)[Particle::IncomingBeam2][0].get().momentum() -
                                                          momenta_.at(1));
}

double MadGraphCollinearProcessBuilder::computeWeight() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphCollinearProcessBuilder:eval") << "Process not properly linked!";

  double weight = 0.;
  momenta_ = mg5_proc_->generateMomenta(sqs_ * std::sqrt(x1_ * x2_), weight);
  mg5_proc_->setMomenta(momenta_);

  CG_DEBUG_LOOP("MadGraphCollinearProcessBuilder:eval")
      << "Particles content:\n"
      << "incoming: x1/2=" << x1_ << "/" << x2_ << ", p4=" << momenta_.at(0) << "/" << momenta_.at(1) << ".";

  return (*flux_)(x1_) * (*flux_)(x2_)*mg5_proc_->eval();
}

ParametersDescription MadGraphCollinearProcessBuilder::description() {
  auto desc = proc::Process::description();
  desc += MadGraphProcessBuilder::description();
  desc.setDescription("MadGraph_aMC process builder (collinear fluxes)");
  desc.add<ParametersDescription>("collinearFluxes", ParametersDescription().setName<std::string>("BudnevEPAProton"));
  return desc;
}

REGISTER_PROCESS("mg5_aMC_coll", MadGraphCollinearProcessBuilder)
