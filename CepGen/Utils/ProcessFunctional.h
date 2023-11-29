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

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Timer.h"

namespace cepgen {
  namespace proc {
    class ProcessFunctional : public SteeredObject<ProcessFunctional> {
    public:
      explicit ProcessFunctional(const ParametersList&, FactorisedProcess&);

      static ParametersDescription description();

      double operator()(double x1, double x2);

    private:
      const std::unique_ptr<Integrator> integr_;
      const std::unique_ptr<ProcessIntegrand> integrand_;
    };
  }  // namespace proc
}  // namespace cepgen
