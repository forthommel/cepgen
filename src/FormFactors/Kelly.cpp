/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen::formfac {
  /// \cite Kelly:2004hm
  class Kelly : public Parameterisation {
  public:
    explicit Kelly(const ParametersList& params)
        : Parameterisation(params),
          ae_(steer<std::vector<double> >("aE")),
          be_(steer<std::vector<double> >("bE")),
          am_(steer<std::vector<double> >("aM")),
          bm_(steer<std::vector<double> >("bM")) {}

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Kelly");
      desc.add("aE", std::vector{1., -0.24});
      desc.add("bE", std::vector{10.98, 12.82, 21.97});
      desc.add("aM", std::vector{1., 0.12});
      desc.add("bM", std::vector{10.97, 18.86, 6.55});
      return desc;
    }

  protected:
    static double computeFF(double tau, const std::vector<double>& as, const std::vector<double>& bs) {
      double num{0.}, den{1.};
      for (size_t i = 0; i < as.size(); ++i)
        num += as.at(i) * std::pow(tau, i);
      for (size_t i = 0; i < bs.size(); ++i)
        den += bs.at(i) * std::pow(tau, i + 1);
      return num / den;
    }

  private:
    void eval() override {
      const auto tau_value = tau(q2_);
      setGEGM(computeFF(tau_value, ae_, be_), MU * computeFF(tau_value, am_, bm_));
    }

    const std::vector<double> ae_, be_;

  protected:
    const std::vector<double> am_, bm_;
  };

  // Alberico et al. fits
  struct Alberico1 final : Kelly {
    using Kelly::Kelly;
    static ParametersDescription description() {
      auto desc = Kelly::description();
      desc.setDescription("Alberico (fit #1)");
      desc.add("aE", std::vector<double>{});
      desc.add("bE", std::vector<double>{});
      desc.add("aM", std::vector{1., 1.53});
      desc.add("bM", std::vector{12.87, 29.16, 41.4});
      return desc;
    }
    void eval() override {
      static constexpr double kC0 = 1.02, kC1 = -0.13;
      const auto ff_m = computeFF(tau(q2_), am_, bm_);
      setGEGM((kC0 + kC1 * q2_) * ff_m, MU * ff_m);
    }
  };
  struct Alberico2 final : Kelly {
    using Kelly::Kelly;
    static ParametersDescription description() {
      auto desc = Kelly::description();
      desc.setDescription("Alberico (fit #2)");
      desc.add("aE", std::vector{1., -0.19});
      desc.add("bE", std::vector{11.12, 15.16, 21.25});
      desc.add("aM", std::vector{1., 1.09});
      desc.add("bM", std::vector{12.31, 25.57, 30.61});
      return desc;
    }
  };
}  // namespace cepgen::formfac
using cepgen::formfac::Kelly;
REGISTER_FORMFACTORS("Kelly", Kelly);
using cepgen::formfac::Alberico1;
using cepgen::formfac::Alberico2;
REGISTER_FORMFACTORS("Alberico1", Alberico1);
REGISTER_FORMFACTORS("Alberico2", Alberico2);
