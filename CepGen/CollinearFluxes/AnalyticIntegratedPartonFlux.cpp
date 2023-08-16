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
#include <memory>

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  class AnalyticIntegratedPartonFlux final : public IntegratedPartonFlux {
  public:
    explicit AnalyticIntegratedPartonFlux(const ParametersList& params)
        : IntegratedPartonFlux(params),
          integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
          unint_flux_(PartonFluxFactory::get().buildCollinearFlux(steer<ParametersList>("unintegratedFlux"))),
          log_q2_(steer<bool>("logQ2")),
          fast_expo_(steer<bool>("fastExpo")) {
      CG_DEBUG("AnalyticIntegratedPartonFlux")
          << "Analytic integrated parton flux evaluator built with:\n"
          << " * integrator: " << integr_->parameters() << "\n"
          << " * unintegrated flux: " << unint_flux_->parameters() << "\n"
          << " * Q^2 virtuality range: " << q2_range_ << " GeV^2 (logarithmic integration? " << log_q2_ << ")";
    }

    static ParametersDescription description() {
      auto desc = IntegratedPartonFlux::description();
      desc.setDescription("Integr.parton flux");
      desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"))
          .setDescription("numerical integrator algorithm to use");
      desc.add<ParametersDescription>("unintegratedFlux", ParametersDescription().setName<std::string>("coll.EPAFlux"))
          .setDescription("(x, Q^2)-dependent flux to integrate with respect to Q^2");
      desc.add<Limits>("q2range", {0., 20.});
      desc.add<bool>("logQ2", true)
          .setDescription("integrate vs. log(Q^2) instead of Q^2 (useful for strongly peaking spectra)");
      desc.add<bool>("fastExpo", false).setDescription("use a fast version of the exponential numerical evaluation");
      return desc;
    }

    double flux(double x) const override {
      Limits q2range;
      if (!computeQ2range(x, q2range))
        return 0.;
      if (log_q2_)
        return integr_->integrate(
            [this, &x](double logq2) {
              const auto q2 = fast_expo_ ? expo(logq2) : std::exp(logq2);
              return unint_flux_->fluxQ2(x, q2) * q2;
            },
            Limits{std::log(q2range.min()), std::log(q2range.max())});
      return integr_->integrate([this, &x](double q2) { return unint_flux_->fluxQ2(x, q2); }, q2range);
    }

    bool fragmenting() const override { return unint_flux_->fragmenting(); }
    pdgid_t partonPdgId() const override { return unint_flux_->partonPdgId(); }
    double mass2() const override { return unint_flux_->mass2(); }

  protected:
    static float expo(float x) {  // cubic spline approximation
      union {
        float f;
        int32_t i;
      } reinterpreter;
      reinterpreter.i = (int32_t)(12102203.0f * x) + 127 * (1 << 23);
      int32_t m = (reinterpreter.i >> 7) & 0xFFFF;  // copy mantissa
      // empirical values for small maximum relative error (8.34e-5):
      reinterpreter.i += ((((((((1277 * m) >> 14) + 14825) * m) >> 14) - 79749) * m) >> 11) - 626;
      return reinterpreter.f;
    }
    const std::unique_ptr<AnalyticIntegrator> integr_;
    const std::unique_ptr<CollinearFlux> unint_flux_;
    const bool log_q2_, fast_expo_;
  };
}  // namespace cepgen
REGISTER_FLUX("integ.IntegratedFlux", AnalyticIntegratedPartonFlux);
