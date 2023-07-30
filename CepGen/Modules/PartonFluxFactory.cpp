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

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"

namespace cepgen {
  std::unique_ptr<KTFlux> PartonFluxFactory::buildKTFlux(const ParametersList& params) const {
    return std::unique_ptr<KTFlux>(dynamic_cast<KTFlux*>(BasePartonFluxFactory::build(params).release()));
  }

  std::unique_ptr<KTFlux> PartonFluxFactory::buildKTFlux(const std::string& name, const ParametersList& params) const {
    return std::unique_ptr<KTFlux>(dynamic_cast<KTFlux*>(BasePartonFluxFactory::build(name, params).release()));
  }

  std::unique_ptr<CollinearFlux> PartonFluxFactory::buildCollinearFlux(const ParametersList& params) const {
    return std::unique_ptr<CollinearFlux>(dynamic_cast<CollinearFlux*>(BasePartonFluxFactory::build(params).release()));
  }

  std::unique_ptr<CollinearFlux> PartonFluxFactory::buildCollinearFlux(const std::string& name,
                                                                       const ParametersList& params) const {
    return std::unique_ptr<CollinearFlux>(
        dynamic_cast<CollinearFlux*>(BasePartonFluxFactory::build(name, params).release()));
  }

  std::unique_ptr<IntegratedPartonFlux> PartonFluxFactory::buildIntegratedFlux(const ParametersList& params) const {
    return std::unique_ptr<IntegratedPartonFlux>(
        dynamic_cast<IntegratedPartonFlux*>(BasePartonFluxFactory::build(params).release()));
  }

  std::unique_ptr<IntegratedPartonFlux> PartonFluxFactory::buildIntegratedFlux(const std::string& name,
                                                                               const ParametersList& params) const {
    return std::unique_ptr<IntegratedPartonFlux>(
        dynamic_cast<IntegratedPartonFlux*>(BasePartonFluxFactory::build(name, params).release()));
  }
}  // namespace cepgen
