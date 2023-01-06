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

#include <Pythia8/Pythia.h>
#include <Pythia8/SigmaEW.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"

namespace cepgen {
  class PythiaProcessBuilder : public proc::Process {
  public:
    explicit PythiaProcessBuilder(const ParametersList& params)
        : Process(params), pythia_(new Pythia8::Pythia), proc_name_(steer<std::string>("processName")) {
      if (proc_name_ == "gmgm2ffbar") {
        const auto pair = steer<ParticleProperties>("pair").pdgid;
        proc_.reset(new Pythia8::Sigma2gmgm2ffbar((int)pair, 0));
      } else
        throw CG_FATAL("PythiaProcessBuilder") << "Invalid process name specified: '" << proc_name_ << "'.";
      proc_->initProc();
    }

    proc::ProcessPtr clone() const {
      throw CG_FATAL("PythiaProcessBuilder") << "Pythia 8 process builder currently does not support cloning.";
    }

    static ParametersDescription description() {
      auto desc = Process::description();
      desc.setDescription("Pythia 8 process handler");
      desc.add<std::string>("processName", "gmgm2ll").setDescription("name of the Pythia 8 process to evaluate");
      return desc;
    }

    double computeWeight() override {
      proc_->sigmaKin();
      return proc_->sigmaHat();
    }

    void addEventContent() override {}

    void fillKinematics() override {}

  private:
    const std::unique_ptr<Pythia8::Pythia> pythia_;
    const std::string proc_name_;

    std::unique_ptr<Pythia8::SigmaProcess> proc_;
  };
}  // namespace cepgen

REGISTER_PROCESS("pythia8", PythiaProcessBuilder);
