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

#include <cmath>

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  class DreesZeppenfeldIntegratedPartonFlux : public IntegratedPartonFlux {
  public:
    explicit DreesZeppenfeldIntegratedPartonFlux(const ParametersList& params) : IntegratedPartonFlux(params) {}

    double mass2() const override { return mp2_; }
    bool fragmenting() const override final { return true; }
    pdgid_t partonPdgId() const override final { return PDG::photon; }

    static ParametersDescription description() {
      auto desc = IntegratedPartonFlux::description();
      desc.setDescription("Drees/Zeppenfeld integ.flux");
      desc.add<pdgid_t>("pdgId", PDG::proton).setDescription("nucleon PDG id");
      return desc;
    }

    double flux(double x) const override {
      Limits q2range;
      if (!computeQ2range(x, q2range))
        return 0.;
      const auto a = 1. + 0.71 / q2range.min(), inv_a = 1. / a;
      return prefactor_ * (1. - x + 0.5 * x * x) *
             (std::log(a) - 11. / 6 + 3. * inv_a - 1.5 * inv_a * inv_a + inv_a * inv_a * inv_a / 3.) / x;
    }
  };
}  // namespace cepgen

REGISTER_INTEGRATED_PARTON_FLUX("DZ", DreesZeppenfeldIntegratedPartonFlux);
