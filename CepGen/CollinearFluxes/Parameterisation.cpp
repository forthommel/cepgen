/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace collflux {
    Parameterisation::Parameterisation(const ParametersList& params)
        : PartonFlux(params), q2_range_(steer<Limits>("q2range")), qscale_(steer<double>("qscale")) {
      q2_range_.min() = std::max(q2_range_.min(), 0.);
    }

    ParametersDescription Parameterisation::description() {
      auto desc = PartonFlux::description();
      desc.setDescription("Unnamed collinear flux");
      desc.add<Limits>("q2range", {0., 1.e4});
      desc.add<double>("qscale", 0.71);
      return desc;
    }

    double Parameterisation::operator()(double /*x*/, double /*kt2*/, double /*mf2*/) const {
      throw CG_FATAL("Parameterisation") << "Calling kt-dependent operator on collinear flux.";
      //return (*this)(x, mf2);
    }
  }  // namespace collflux
}  // namespace cepgen