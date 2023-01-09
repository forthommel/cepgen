/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process2to4.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcessBuilder.h"

using namespace cepgen;

class MadGraphKTProcessBuilder : public MadGraphProcessBuilder, public proc::Process2to4 {
public:
  explicit MadGraphKTProcessBuilder(const ParametersList&);
  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphKTProcessBuilder(*this)); }

  static ParametersDescription description();

  void prepareProcessKinematics() override;
  double computeCentralMatrixElement() const override;
};

MadGraphKTProcessBuilder::MadGraphKTProcessBuilder(const ParametersList& params)
    : MadGraphProcessBuilder(params), Process2to4(params, std::array<pdgid_t, 2>{}, 0) {
  const auto& interm_part = mg5_proc_->intermediatePartons();
  const auto& cent_sys = mg5_proc_->centralSystem();
  CG_DEBUG("MadGraphKTProcessBuilder") << "MadGraph_aMC process created for:\n\t"
                                       << "* interm. parts.: " << interm_part << "\n\t"
                                       << "* central system: " << cent_sys << ".";
  setIntermediatePartons({(pdgid_t)interm_part[0], (pdgid_t)interm_part[1]});
  setProducedParticles(std::vector<pdgid_t>(cent_sys.begin(), cent_sys.end()));
}

void MadGraphKTProcessBuilder::prepareProcessKinematics() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphKTProcessBuilder") << "Process not properly linked!";

  CG_INFO("MadGraphKTProcessBuilder") << "Preparing process kinematics for card at \"" << params_card_ << "\".";
  mg5_proc_->initialise(params_card_);
}

double MadGraphKTProcessBuilder::computeCentralMatrixElement() const {
  if (!mg5_proc_)
    CG_FATAL("MadGraphKTProcessBuilder:eval") << "Process not properly linked!";

  CG_DEBUG_LOOP("MadGraphKTProcessBuilder:eval")
      << "Particles content:\n"
      << "incoming: " << q1_ << " (m=" << q1_.mass() << "), " << q2_ << " (m=" << q2_.mass() << ")\n"
      << "outgoing: " << p_c1_ << " (m=" << p_c1_.mass() << "), " << p_c2_ << " (m=" << p_c2_.mass() << ").";
  mg5_proc_->setMomentum(0, q1_);    // first incoming parton
  mg5_proc_->setMomentum(1, q2_);    // second incoming parton
  mg5_proc_->setMomentum(2, p_c1_);  // first outgoing central particle
  mg5_proc_->setMomentum(3, p_c2_);  // second outgoing central particle

  return mg5_proc_->eval();
}

ParametersDescription MadGraphKTProcessBuilder::description() {
  auto desc = Process2to4::description();
  desc += MadGraphProcessBuilder::description();
  desc.setDescription("MadGraph_aMC process builder (kT-factorised ME)");
  return desc;
}

REGISTER_PROCESS("mg5_aMC", MadGraphKTProcessBuilder)
