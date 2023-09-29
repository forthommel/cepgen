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
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// Regge-like NNLO, low-x F2/FL structure functions
    /// \cite Baruah:2014mwa
    class BaruahDasSarma final : public Parameterisation {
    public:
      explicit BaruahDasSarma(const ParametersList& params)
          : Parameterisation(params),
            embedded_gluon_evol_(steer<bool>("embeddedGluonEvolution")),
            integr_g_(embedded_gluon_evol_ ? AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))
                                           : nullptr),
            integr_f2_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
            alpha_s_(AlphaSFactory::get().build(steer<ParametersList>("alphaS"))),
            gluon_flux_(embedded_gluon_evol_
                            ? nullptr
                            : CollinearFluxFactory::get().build(steer<ParametersList>("externalGluonFlux"))),
            nf_(steer<int>("numFlavours")),
            esq_(steer<double>("eSquare")),
            lambda_g_(steer<double>("lambdaG")),
            beta0_(11. - 2. * nf_ / 3.) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Baruah-Das-Sarma");
        desc.add<bool>("embeddedGluonEvolution", false);
        desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"))
            .setDescription("Steering parameters for the analytical integrator");
        desc.add<ParametersDescription>("alphaS", AlphaSFactory::get().describeParameters("pegasus"))
            .setDescription("strong coupling evolution algorithm");
        desc.add<ParametersDescription>("externalGluonFlux",
                                        ParametersDescription(ParametersList()
                                                                  .setName<std::string>("LHAPDF")
                                                                  .set<std::string>("set", "MSTW2008nnlo90cl")
                                                                  .set<pdgid_t>("partonPdgId", PDG::gluon)));
        desc.add<int>("numFlavours", 3).setDescription("number of active light flavours");
        desc.add<double>("eSquare", 5. / 18).setDescription("average squared charge for even number of flavours");
        desc.add<double>("lambdaG", 1.08).setDescription("Regge gluon distribution function exponent");
        return desc;
      }

      void eval() override {
        const auto delta = [](double x) { return x == 0.; };
        auto c2_1 = [&](double w) {
          const auto w1 = 1. - w, l0 = std::log(w), l1 = std::log(w1);
          return nf_ * ((2. - 4. * w * w1) * (l1 - l0) - 2. + 16. * w * w1);
        };
        auto c2_2 = [&](double w) {
          const auto w1 = 1. - w, l0 = std::log(w), l1 = std::log(w1);
          return nf_ * ((6.445 + 209.4 * (1. - w)) * l1 * l1 * l1 - 24. * l1 * l1 + (149. / w - 1483.) * l1 +
                        l1 * l0 * (-871.8 * l1 - 724.1 * l0) + 5.319 * l0 * l0 * l0 - 59.48 * l0 * l0 - 284.8 * l0 +
                        11.9 / w + 392.4 - 0.28 * delta(1. - w));
        };
        auto c2_3 = [&](double w) {
          const auto w1 = 1. - w, l0 = std::log(w), l1 = std::log(w1), fl11 = 1. /*FIXME charge factor*/;
          const auto x = args_.xbj, x1 = x /*FIXME*/;
          return nf_ * (966. / 81 * l1 * l1 * l1 * l1 * l1 - 1871. / 18 * l1 * l1 * l1 * l1 + 89.31 * l1 * l1 * l1 +
                        979.2 * l1 * l1 - 2405. * l1 + 1372. * x1 * l1 * l1 * l1 * l1 - 15729. - 310'510. * x +
                        331'570 * x * x - 244'150 * x * l0 * l0 - 253.3 * x * l0 * l0 * l0 * l0 * l0 +
                        l0 * l1 * (138'230. - 237'010 * l0) - 11'860 * l0 - 700.8 * l0 * l0 - 1440 * l0 * l0 * l0 +
                        4951. / 162 * l0 * l0 * l0 * l0 - 134. / 9 * l0 * l0 * l0 * l0 * l0 -
                        (6362.54 - 932.089 * l0) / x + 0.625 * delta(x1)) +
                 nf_ * nf_ *
                     (131. / 81 * l1 * l1 * l1 * l1 - 14.72 * l1 * l1 * l1 + 3.607 * l1 * l1 - 226.1 * l1 + 4.762 -
                      190 * x - 818.4 * x * x - 4019. * x * l0 * l0 - l0 * l1 * (791.5 + 4646 * l0) + 739. * l0 +
                      418. * l0 * l0 + 104.3 * l0 * l0 * l0 + 809. / 81 * l0 * l0 * l0 * l0 +
                      12. / 9 * l0 * l0 * l0 * l0 * l0 + 84.423 / x) +
                 fl11 * nf_ * nf_ *
                     (3.211 * l1 * l1 + 19.04 * x * l1 + 0.623 * x1 * l1 * l1 * l1 - 64.47 * x + 121.6 * x * x -
                      45.82 * x * x * x - x * l0 * l1 * (31.68 + 37.24 * l0) + 11.27 * x * x * l0 * l0 * l0 -
                      82.4 * x * l0 - 16.08 * x * l0 * l0 + 520. / 81 * x * l0 * l0 * l0 +
                      20. / 27 * x * l0 * l0 * l0 * l0);
        };
        auto c2 = [this, &c2_1, &c2_2, &c2_3](double w, double q2) {
          const auto prefactor = 0.25 * (*alpha_s_)(std::sqrt(q2)) * M_1_PI;
          return prefactor * c2_1(w) + prefactor * prefactor * c2_2(w) + prefactor * prefactor * prefactor * c2_3(w);
        };
        setF2(args_.xbj * esq_ *
              integr_f2_->integrate(
                  [this, &c2](double w) -> double {
                    return c2(w, args_.q2) * (embedded_gluon_evol_
                                                  ? nnloGluonDistribution(args_.xbj / w, args_.q2)
                                                  : (gluon_flux_->fluxQ2(args_.xbj / w, args_.q2) /*/ args_.xbj*/));
                  },
                  Limits{args_.xbj, 1.}));
      }

    private:
      double nnloGluonDistribution(double x, double q2) const {
        auto p = [this](double x) {
          return 12. / beta0_ *
                 ((11. / 12 - nf_ * 1. / 18) + std::log(1. - x) +
                  integr_g_->integrate(
                      [this](double w) {
                        return (std::pow(w, 1. + lambda_g_) - 1.) / (1. - w) +
                               std::pow(w, lambda_g_) * (w * (1. - w) + (1. - w) / w);
                      },
                      Limits{x, 1.}));
        };
        return g_funct_c_ * std::pow(q2, p(x)) / x;
      }
      double g_funct_c_{1.};

      const bool embedded_gluon_evol_;
      const std::unique_ptr<AnalyticIntegrator> integr_g_, integr_f2_;
      const std::unique_ptr<Coupling> alpha_s_;
      const std::unique_ptr<CollinearFlux> gluon_flux_;
      const int nf_;
      const double esq_;
      const double lambda_g_;
      const double beta0_;
    };
  }  // namespace strfun
}  // namespace cepgen
using cepgen::strfun::BaruahDasSarma;
REGISTER_STRFUN(107, BaruahDasSarma);
