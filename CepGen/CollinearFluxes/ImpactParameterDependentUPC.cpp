/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include <gsl/gsl_sf_bessel.h>

#include <cmath>

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/FunctionsWrappers.h"
#include "CepGen/Utils/Limits.h"

//FIXME
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"

namespace cepgen {
  class ImpactParameterDependentUPC : public IntegratedPartonFlux {
  public:
    explicit ImpactParameterDependentUPC(const ParametersList& params)
        : IntegratedPartonFlux(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          q2_range_(steer<Limits>("q2range")),
          integr_ff_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("analyticalIntegrator"))),
          integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("analyticalIntegrator"))),
          flux_args_(new FluxArguments{0., 0.}),
          mass2_(hi_.mass() * hi_.mass()),
          prefactor_((int)hi_.Z * (int)hi_.Z * constants::ALPHA_EM * M_1_PI * M_1_PI) {
      static const auto rho = [this](double r) {
        static const double rA = 6.63, delta = 0.459;  // fm
        return rho0_ / (1. + exp((r - rA) / delta));
      };
      const Limits r_range{0., 30.};  // r integration range (in fm)

      {  // extract normalisation for the r range considered
        rho0_ = 1. / integr_ff_->integrate([&](double r) { return rho(r); }, r_range);
        const auto norm = integr_ff_->integrate([&](double r) { return rho(r); }, r_range);
        if (fabs(norm - 1.) > 1.e-12)
          throw CG_FATAL("ImpactParameterDependentUPC")
              << "Failed to normalise the Wood-Saxon charge distribution rho(r)! Expecting 1, got 1+-" << (norm - 1.)
              << ".";
      }
      auto nucl_form_fac = [&](double q2) -> double {
        if (q2 == 0.)
          return 0.;
        const auto q = std::sqrt(q2);
        return 4. * M_PI * integr_ff_->integrate([&q](double r) { return r * rho(r) * sin(q * r) / q; }, r_range);
      };
      {
        utils::Graph1D gr_ff;
        for (const auto& x : q2_range_.generate(1000))
          gr_ff.addPoint(x, nucl_form_fac(x));
        DrawerFactory::get()
            .build("root", ParametersList().set<std::string>("filename", "nucl_form_factor"))
            ->draw(gr_ff, utils::Drawer::Mode::grid);
      }

      // initialise the function to integrate
      func_.reset(new utils::Function1D([&](double kt, void* params) {
        if (!params)
          throw CG_FATAL("ImpactParameterDependentUPC") << "Invalid parameters block fed to the integrand!";
        const auto& args = *static_cast<FluxArguments*>(params);
        const double kt2 = kt * kt, m2 = kt2 + args.ene2;
        return kt2 * nucl_form_fac(m2) / m2 * gsl_sf_bessel_J1(args.b * kt);
      }));
      auto pointlike_form_fac = [&](double b) -> double {
        if (b == 0.)
          return 0.;
        const double xi = 1.e-3;
        const double gamma = 2676.;
        const double omega = xi * 5020;
        const auto k = b * omega / gamma;
        return prefactor_ * k * k / omega / b / b *
               (std::pow(gsl_sf_bessel_K1(k), 2) + std::pow(gsl_sf_bessel_K0(k) / gamma, 2));
      };

      {
        utils::Graph1D gr_ff, gr_ff_pl;
        for (const auto& b : Limits(0., 30.).generate(500)) {
          FluxArguments args;
          args.b = b;
          args.ene2 = 14.e3 * 14.e3;
          //gr_ff.addPoint(b, integr_->integrate(*func_, args, q2_range_));
          CG_LOG << b << ":" << pointlike_form_fac(b);
          gr_ff_pl.addPoint(b, pointlike_form_fac(b));
        }
        DrawerFactory::get().build("root")->draw(
            {&gr_ff_pl}, "photon_number", "", utils::Drawer::Mode::logy | utils::Drawer::Mode::grid);
      }
      //exit(0);
      CG_INFO("ImpactParameterDependentUPC")
          << "IP-dependent UPC collinear flux evaluator initialised.\n\t"
          << "Q^2 integration range: " << q2_range_ << " GeV^2\n\t"
          << "r integration range: " << r_range << " fm, rho0 normalisation: " << rho0_ << "\n\t"
          << "Nucleon/HI: " << hi_ << ".";
    }

    bool fragmenting() const override final { return false; }
    pdgid_t partonPdgId() const override final { return PDG::photon; }

    static ParametersDescription description() {
      auto desc = IntegratedPartonFlux::description();
      desc.setDescription("Impact param.-dependent UPC photon flux");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::proton());
      desc.add<Limits>("q2range", {0., 1.e4}).setDescription("kinematic range for the parton virtuality, in GeV^2");
      desc.add<ParametersDescription>("analyticalIntegrator", ParametersDescription().setName<std::string>("gsl"))
          .setDescription("Steering parameters for the analytical integrator");
      return desc;
    }

    double mass2() const override { return mass2_; }

    double flux(double x) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      flux_args_->b = 0.5;
      flux_args_->ene2 = 0.25 * x * x * 14.e3 * 14.e3;
      const auto q2_integr_flux = integr_->integrate(*func_, (void*)flux_args_.get(), q2_range_);
      CG_LOG << x << ":" << q2_integr_flux;
      return prefactor_ * q2_integr_flux * q2_integr_flux;
    }

  private:
    const HeavyIon hi_;
    const Limits q2_range_;
    std::unique_ptr<utils::Function1D> func_;
    std::unique_ptr<AnalyticIntegrator> integr_ff_;
    std::unique_ptr<AnalyticIntegrator> integr_;
    struct FluxArguments {
      double b{0.}, ene2{0.};
    };
    std::unique_ptr<FluxArguments> flux_args_;
    const double mass2_;
    const double prefactor_{0.};
    double rho0_{1.};
  };
}  // namespace cepgen
REGISTER_INTEGRATED_PARTON_FLUX("ImpactParameterDependentUPC", ImpactParameterDependentUPC);
