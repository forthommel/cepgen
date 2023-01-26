/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2023  Laurent Forthomme
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

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/PartonFlux.h"

namespace cepgen {
  class EPACollinearFlux : public PartonFlux {
  public:
    explicit EPACollinearFlux(const ParametersList& params)
        : PartonFlux(params), yrange_(steer<Limits>("yRange")), dyrange_(steer<Limits>("dyRange")) {}

    double operator()(double x, double mf2) const override { return 0.; }

    double operator()(double x, double kt2, double mf2) const override {
      (void)kt2;
      return operator()(x, mf2);
    }

    static ParametersDescription description() {
      auto desc = PartonFlux::description();
      desc.add<Limits>("q2Range", {0., 1.e6});
      desc.add<Limits>("wRange", {0., 100.});
      desc.add<Limits>("yRange", {0., 1.});
      desc.add<Limits>("dyRange", {0., 1.});
      return desc;
    }

  private:
    const Limits yrange_;
    const Limits dyrange_;
  };
}  // namespace cepgen

REGISTER_FLUX("EPACollinear", EPACollinearFlux);
