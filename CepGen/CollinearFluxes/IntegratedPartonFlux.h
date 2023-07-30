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

#ifndef CepGen_CollinearFluxes_IntegratedPartonFlux_h
#define CepGen_CollinearFluxes_IntegratedPartonFlux_h

#include <memory>

#include "CepGen/CollinearFluxes/CollinearFlux.h"

namespace cepgen {
  class AnalyticIntegrator;
  class IntegratedPartonFlux : public PartonFlux {
  public:
    explicit IntegratedPartonFlux(const ParametersList&);

    static ParametersDescription description();

    /// Compute the collinear flux for this x value
    virtual double flux(double x) const = 0;

    bool ktFactorised() const override final { return false; }
    bool integratedQ2() const override final { return true; }

  protected:
    /// Integration range for the flux
    const Limits q2_range_{0., 1.e4};
  };
}  // namespace cepgen

#endif
