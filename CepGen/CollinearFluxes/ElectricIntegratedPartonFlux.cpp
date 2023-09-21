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
  class GenericElectricIntegratedPartonFlux : public IntegratedPartonFlux {
  public:
    explicit GenericElectricIntegratedPartonFlux(const ParametersList& params)
        : IntegratedPartonFlux(params),
          q2_scale_(steer<double>("scale")),
          polyn_coeff_(steer<std::vector<double> >("polynomialCoefficients")),
          a_taming_(steer<bool>("aTaming")) {}

    double mass2() const override { return mp2_; }
    bool fragmenting() const override final { return true; }
    pdgid_t partonPdgId() const override final { return PDG::photon; }

    static ParametersDescription description() {
      auto desc = IntegratedPartonFlux::description();
      desc.setDescription("Generic electric integ.flux");
      desc.add<pdgid_t>("pdgId", PDG::proton).setDescription("nucleon PDG id");
      desc.add<double>("scale", 0.71)
          .setDescription("scaling (in GeV^2) (0.71 for r_p = 0.81 fm, 0.66 for r_p = 0.84 fm)");
      desc.add<std::vector<double> >("polynomialCoefficients", {-11. / 6, 3., -3. / 2, 1. / 3});
      desc.add<bool>("aTaming", false);
      return desc;
    }

    double flux(double x) const override {
      Limits q2range;
      if (!computeQ2range(x, q2range))
        return 0.;
      const auto a = 1. + q2_scale_ / q2range.min(), inv_a = 1. / a;
      double polyn = 0.;
      for (size_t i = 0; i < polyn_coeff_.size(); ++i)
        polyn += polyn_coeff_.at(i) * std::pow(inv_a, i);
      const auto factor_a = a_taming_ ? (a + 3.) / (a - 1.) : 1.;
      return prefactor_ * (1. - x + 0.5 * x * x) / x * (factor_a * std::log(a) + polyn);
    }

  private:
    const double q2_scale_;
    const std::vector<double> polyn_coeff_;
    const bool a_taming_;
  };

  struct DreesZeppenfeldIntegratedPartonFlux final : public GenericElectricIntegratedPartonFlux {
    using GenericElectricIntegratedPartonFlux::GenericElectricIntegratedPartonFlux;
    static ParametersDescription description() {
      auto desc = GenericElectricIntegratedPartonFlux::description();
      desc.setDescription("Drees/Zeppenfeld integ.flux");
      desc.add<std::vector<double> >("polynomialCoefficients", {-11. / 6, 3., -3. / 2, 1. / 3});
      desc.add<bool>("aTaming", false);
      return desc;
    }
  };

  struct ElectricIntegratedPartonFlux final : public GenericElectricIntegratedPartonFlux {
    using GenericElectricIntegratedPartonFlux::GenericElectricIntegratedPartonFlux;
    static ParametersDescription description() {
      auto desc = GenericElectricIntegratedPartonFlux::description();
      desc.setDescription("Electric integ.flux");
      desc.add<std::vector<double> >("polynomialCoefficients", {-17. / 6, -4. / 3, 1. / 6});
      desc.add<bool>("aTaming", true);
      return desc;
    }
  };
}  // namespace cepgen

REGISTER_INTEGRATED_PARTON_FLUX("DZ", DreesZeppenfeldIntegratedPartonFlux);
REGISTER_INTEGRATED_PARTON_FLUX("Electric", ElectricIntegratedPartonFlux);
