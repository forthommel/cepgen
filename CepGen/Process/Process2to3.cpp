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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/PhaseSpaceGenerator.h"
#include "CepGen/Process/Process2to3.h"

namespace cepgen {
  namespace proc {
    Process2to3::Process2to3(const ParametersList& params, std::array<pdgid_t, 2> partons, pdgid_t cs_id)
        : KTProcess(params, partons, {cs_id}), single_limits_(params) {}

    void Process2to3::setCuts(const cuts::Central& single) { single_limits_ = single; }

    void Process2to3::preparePhaseSpace() {
      pgen_->initialise();
      prepareProcessKinematics();
    }

    double Process2to3::computeKTFactorisedMatrixElement() {
      if (!pgen_->generate())
        return 0.;
      // compute the central 2-to-1 matrix element
      const double amat2 = computeCentralMatrixElement();
      if (amat2 <= 0.)  // skip computing the fluxes if no contribution
        return 0.;

      //=================================================================
      // factor 1/4 from jacobian of transformations
      // factors 1/pi and 1/pi due to integration over
      //     d^2(kappa_1)d^2(kappa_2) instead of d(kappa_1^2)d(kappa_2^2)
      //=================================================================

      return amat2 * pow(4. * x1_ * x2_ * s_ * M_PI, -2) * 0.25 * constants::GEVM2_TO_PB * qt1_ * qt2_;
    }

    void Process2to3::fillCentralParticlesKinematics() {
      //--- outgoing central particle
      auto& oc = (*event_)[Particle::CentralSystem][0].get();
      oc.setStatus(Particle::Status::Undecayed);
    }

    //----- utilities

    ParametersDescription Process2to3::description() {
      auto desc = KTProcess::description();
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
