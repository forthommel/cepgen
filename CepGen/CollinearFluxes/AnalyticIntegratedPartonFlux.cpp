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

#include <memory>

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"

namespace cepgen {
  class AnalyticIntegratedPartonFlux final : public IntegratedPartonFlux {
  public:
    explicit AnalyticIntegratedPartonFlux(const ParametersList& params)
        : IntegratedPartonFlux(params),
          integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
          unint_flux_(PartonFluxFactory::get().buildCollinearFlux(steer<ParametersList>("unintegratedFlux"))) {}

    static ParametersDescription description() {
      auto desc = IntegratedPartonFlux::description();
      desc.setDescription("Integrated parton flux evaluator");
      desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"));
      desc.add<ParametersDescription>("unintegratedFlux",
                                      ParametersDescription().setName<std::string>("coll.PlainEPA"));
      return desc;
    }

    double flux(double x) const override {
      return integr_->integrate([this, &x](double q2) { return unint_flux_->fluxQ2(x, q2); }, q2_range_);
    }

    bool fragmenting() const override { return unint_flux_->fragmenting(); }
    pdgid_t partonPdgId() const override { return unint_flux_->partonPdgId(); }
    double mass2() const override { return unint_flux_->mass2(); }

  protected:
    const std::unique_ptr<AnalyticIntegrator> integr_;
    const std::unique_ptr<CollinearFlux> unint_flux_;
  };
}  // namespace cepgen
REGISTER_FLUX("integ.IntegratedFlux", AnalyticIntegratedPartonFlux);
