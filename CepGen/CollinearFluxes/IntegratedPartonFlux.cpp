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

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"

namespace cepgen {
  IntegratedPartonFlux::IntegratedPartonFlux(const ParametersList& params)
      : PartonFlux(params), q2_range_(steer<Limits>("Q2range")) {}

  ParametersDescription IntegratedPartonFlux::description() {
    auto desc = PartonFlux::description();
    desc.setDescription("Integrated parton flux evaluator");
    desc.add<Limits>("Q2range", Limits{0., 20.})
        .setDescription("integration bounds for the parton virtuality, in GeV^2");
    return desc;
  }

  bool IntegratedPartonFlux::computeQ2range(double x, Limits& q2limits) const {
    if (!x_range_.contains(x, true))
      return false;
    q2limits = Limits{mass2() * x * x / (1. - x), q2_range_.max()};
    return q2_range_.contains(q2limits.min());
  }
}  // namespace cepgen
