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

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  class DGLAPGluonFlux : public CollinearFlux {
  public:
    explicit DGLAPGluonFlux(const ParametersList& params)
        : CollinearFlux(params),
          integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
          nf_(steer<int>("numFlavours")),
          lambda_g_(steer<double>("lambdaG")),
          beta0_(11. - 2. * nf_ / 3.) {
      //g_func_c_ = 1. / integr_->integrate()
    }

    double fluxQ2(double x, double q2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      return g_func_c_ * std::pow(q2, pFunction(x)) / x;
    }

    static ParametersDescription description() {
      auto desc = CollinearFlux::description();
      desc.setDescription("DGLAP gluon flux");
      desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"))
          .setDescription("Steering parameters for the analytical integrator");
      desc.add<int>("numFlavours", 3).setDescription("number of active light flavours");
      desc.add<double>("lambdaG", 1.08).setDescription("Regge gluon distribution function exponent");
      return desc;
    }

    pdgid_t partonPdgId() const override final { return PDG::gluon; }
    double mass2() const override final { return mp2_; }

  private:
    double pFunction(double x) const {
      return 12. / beta0_ *
             ((11. / 12 - nf_ * 1. / 18) + std::log(1. - x) +
              integr_->integrate(
                  [this](double w) {
                    return (std::pow(w, 1. + lambda_g_) - 1.) / (1. - w) +
                           std::pow(w, lambda_g_) * (w * (1. - w) + (1. - w) / w);
                  },
                  Limits{x, 1.}));
    }
    const std::unique_ptr<AnalyticIntegrator> integr_;
    const int nf_;
    const double lambda_g_;
    const double beta0_;
    double g_func_c_{1.};
  };
}  // namespace cepgen
REGISTER_COLLINEAR_FLUX("DGLAP", DGLAPGluonFlux);
