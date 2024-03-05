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
#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;

class ChiliPhaseSpaceGenerator : public cepgen::PhaseSpaceGenerator {
public:
  explicit ChiliPhaseSpaceGenerator(const ParametersList& params)
      : PhaseSpaceGenerator(params),
        part_psgen_(PartonsPhaseSpaceGeneratorFactory::get().build(steer<std::string>("partonsGenerator"))),
        central_(steer<std::vector<int> >("ids")),
        model_(new chili::Model([](int i, int j) -> std::vector<int> {
          if (i == j && i == 22)
            return {};
          return {};
        })) {
    for (const auto& pdgid : spdgids_t{1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 11, 13, 15, 21, 22, 23}) {
      const auto pdg_info = PDG::get()(pdgid);
      model_->Mass(pdgid) = pdg_info.mass;
      model_->Width(pdgid) = pdg_info.width;
    }
  }

  static ParametersDescription description() {
    auto desc = PhaseSpaceGenerator::description();
    desc.setDescription("Chili central phase space generator");
    return desc;
  }

  pdgids_t partons() const override {
    CG_ASSERT(part_psgen_);
    return pdgids_t{part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
  }

  void setCentral(const std::vector<int>& central) override { central_ = central; }
  std::vector<int> central() const override { return central_; }

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
  const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
  std::vector<int> central_;
  const std::unique_ptr<chili::Model> model_;
  proc::FactorisedProcess* proc_;  //NOT owning

  double central_weight_{0.};
};
REGISTER_PHASE_SPACE_GENERATOR("chili", ChiliPhaseSpaceGenerator);
