/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"

namespace cepgen::helaconia {
  class Process final : public cepgen::proc::Process {
  public:
    explicit Process(const ParametersList& params) : proc::Process(params) {}

    proc::ProcessPtr clone() const override { return std::make_unique<Process>(*this); }

    void addEventContent() override {
      proc::Process::setEventContent({{Particle::Role::IncomingBeam1, {PDG::proton}},
                                      {Particle::Role::IncomingBeam2, {PDG::proton}},
                                      {Particle::Role::Parton1, {PDG::photon}},
                                      {Particle::Role::Parton2, {PDG::photon}},
                                      {Particle::Role::OutgoingBeam1, {PDG::proton}},
                                      {Particle::Role::OutgoingBeam2, {PDG::proton}},
                                      {Particle::Role::CentralSystem, {}}});
    }

    double computeWeight() override { return 0.; }

    void prepareKinematics() override {}

    void fillKinematics() override {}

    static ParametersDescription description() {
      auto desc = proc::Process::description();
      desc.setDescription("HELAC-Onia process");
      return desc;
    }

  private:
  };
}  // namespace cepgen::helaconia
using HelaconiaProcess = cepgen::helaconia::Process;
REGISTER_PROCESS("helaconia", HelaconiaProcess);
