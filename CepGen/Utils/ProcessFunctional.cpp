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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/ProcessFunctional.h"

namespace cepgen {
  namespace proc {
    ProcessFunctional::ProcessFunctional(const ParametersList& params, FactorisedProcess& proc)
        : SteeredObject(params),
          integr_(IntegratorFactory::get().build(steer<ParametersList>("integrator"))),
          integrand_(new ProcessIntegrand(proc)) {}

    ParametersDescription ProcessFunctional::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Helper to compute process as a functional");
      desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("Vegas"))
          .setDescription("Steering parameters for the multidimensional integrator");
      return desc;
    }

    double ProcessFunctional::operator()(double x1, double x2) {
      const auto x_mean = std::sqrt(x1 * x2);
      integrand_->process().kinematics().cuts().remnants.xi = {x_mean, x_mean};
      integrand_->process().initialise();
      return integr_->integrate(*integrand_);
    }
  }  // namespace proc
}  // namespace cepgen
