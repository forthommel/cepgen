/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <Chili/Model/Model.hh>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsCollinearPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsKTPhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;

template <typename T>
class ChiliPhaseSpaceGenerator : public cepgen::PhaseSpaceGenerator {
public:
  explicit ChiliPhaseSpaceGenerator(const ParametersList& params)
      : PhaseSpaceGenerator(params), part_psgen_(new T(params)), central_(steer<std::vector<int> >("ids")) {}

  static ParametersDescription description() {
    auto desc = PhaseSpaceGenerator::description();
    desc.setDescription("Chili central phase space generator");
    return desc;
  }

  pdgids_t partons() const override {
    CG_ASSERT(part_psgen_);
    return pdgids_t{part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
  }

  pdgids_t central() const override { return pdgids_t(central_.begin(), central_.end()); }

  void initialise(proc::FactorisedProcess* process) override {
    CG_ASSERT(part_psgen_);
    part_psgen_->initialise(process);
    proc_ = process;
  }

  bool generate() override {
    CG_ASSERT(part_psgen_);
    if (!part_psgen_->generatePartonKinematics())
      return false;
    return true;
  }

  double weight() const override {
    const auto fluxes_weight = part_psgen_->fluxes();
    if (!utils::positive(fluxes_weight))
      return 0.;
    return fluxes_weight * central_weight_;
  }

  //FIXME only works for 2-to-4
  double that() const override {
    return 0.5 * ((proc_->q1() - proc_->pc(0)).mass2() + (proc_->q2() - proc_->pc(1)).mass2());
  }

  // FIXME only works for 2-to-4
  double uhat() const override {
    return 0.5 * ((proc_->q1() - proc_->pc(1)).mass2() + (proc_->q2() - proc_->pc(0)).mass2());
  }

private:
  std::unique_ptr<chili::Model> model_;
  proc::FactorisedProcess* proc_;  //NOT owning

  const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
  const std::vector<int> central_;

  double central_weight_{0.};
};
using KTChili = ChiliPhaseSpaceGenerator<cepgen::PartonsKTPhaseSpaceGenerator>;
using CollChili = ChiliPhaseSpaceGenerator<cepgen::PartonsCollinearPhaseSpaceGenerator>;
REGISTER_PSGEN("ktchili", KTChili);
REGISTER_PSGEN("collchili", CollChili);
