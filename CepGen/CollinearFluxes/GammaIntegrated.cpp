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

#include <cmath>

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/FunctionsWrappers.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class GammaIntegrated : public CollinearFlux {
  public:
    explicit GammaIntegrated(const ParametersList& params)
        : CollinearFlux(params),
          kt2_range_(steer<Limits>("kt2range")),
          flux_(
              *dynamic_cast<KTFlux*>(PartonFluxFactory::get().build(steer<ParametersList>("ktPartonFlux")).release())),
          integr_(AnalyticIntegratorFactory::get().build(params.get<ParametersList>("analyticalIntegrator"))) {
      if (!flux_.ktFactorised())
        throw CG_FATAL("GammaIntegrated") << "Input flux has to be unintegrated.";
      // initialise the functions to integrate
      func_q2_.reset(new utils::Function1D([&](double kt2, void* params) {
        const auto& args = *static_cast<FluxArguments*>(params);
        return flux_.fluxQ2(args.x, kt2, args.var2);
      }));
      func_mx2_.reset(new utils::Function1D([&](double kt2, void* params) {
        const auto& args = *static_cast<FluxArguments*>(params);
        return flux_.fluxMX2(args.x, kt2, args.var2);
      }));
      CG_INFO("GammaIntegrated") << "kt flux-integrated collinear flux evaluator initialised.\n\t"
                                 << "Q^2 integration range: " << kt2_range_ << " GeV^2\n\t"
                                 << "Unintegrated flux: " << flux_.name() << ".";
    }

    bool fragmenting() const override final { return flux_.fragmenting(); }
    pdgid_t partonPdgId() const override final { return flux_.partonPdgId(); }
    double mass2() const override final { return flux_.mass2(); }

    static ParametersDescription description() {
      auto desc = CollinearFlux::description();
      desc.setDescription("kt-integrated photon flux");
      desc.add<Limits>("q2range", {0., 1.e4}).setDescription("kinematic range for the parton virtuality, in GeV^2");
      desc.add<ParametersDescription>("ktFlux", ParametersDescription().setName<std::string>("BudnevElasticKT"))
          .setDescription("Type of unintegrated kT-dependent parton flux");
      desc.add<ParametersDescription>("analyticalIntegrator", ParametersDescription().setName<std::string>("gsl"))
          .setDescription("Steering parameters for the analytical integrator");
      return desc;
    }

    double fluxQ2(double x, double q2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      return integr_->integrate(*func_q2_, FluxArguments{x, q2}, kt2_range_) / x;
    }

    double fluxMX2(double x, double mx2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      return integr_->integrate(*func_mx2_, FluxArguments{x, mx2}, kt2_range_) / x;
    }

  private:
    const Limits kt2_range_;
    const KTFlux& flux_;
    std::unique_ptr<utils::Function1D> func_q2_, func_mx2_;
    std::unique_ptr<AnalyticIntegrator> integr_;

    struct FluxArguments {
      double x{0.}, var2{0.};
    };
  };
}  // namespace cepgen
REGISTER_FLUX("GammaIntegrated", GammaIntegrated);
