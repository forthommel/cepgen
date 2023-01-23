/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Process/Process2to4.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace proc {
    const Limits Process2to4::x_limits_{0., 1.};

    Process2to4::Process2to4(const ParametersList& params, std::array<pdgid_t, 2> partons, pdgid_t cs_id)
        : KTProcess(params, partons, {cs_id, cs_id}), cs_prop_(PDG::get()(cs_id)), single_limits_(params) {}

    void Process2to4::setCuts(const cuts::Central& single) { single_limits_ = single; }

    void Process2to4::preparePhaseSpace() {
      if (cs_prop_.pdgid == PDG::invalid)  // ensure the central particles properties are correctly initialised
        cs_prop_ = PDG::get()(steer<ParticleProperties>("pair").pdgid);

      pgen_->initialise();

      prepareProcessKinematics();
    }

    double Process2to4::computeKTFactorisedMatrixElement() {
      if (!pgen_->generate())
        return 0.;
      //--- compute the central 2-to-2 matrix element
      const double amat2 = computeCentralMatrixElement();
      if (amat2 <= 0.)  // skip computing the fluxes if no contribution
        return 0.;

      //=================================================================
      // factor 1/4 from jacobian of transformations
      // factors 1/pi and 1/pi due to integration over
      //     d^2(kappa_1)d^2(kappa_2) instead of d(kappa_1^2)d(kappa_2^2)
      //=================================================================

      return amat2 * pow(4. * x1_ * x2_ * s_ * M_PI, -2) * 0.25 * constants::GEVM2_TO_PB * pt_diff_ * qt1_ * qt2_;
    }

    void Process2to4::fillCentralParticlesKinematics() {
      //--- randomise the charge of outgoing system
      short sign = (drand() > 0.5) ? +1 : -1;

      //--- first outgoing central particle
      auto& oc1 = (*event_)[Particle::CentralSystem][0].get();
      oc1.setChargeSign(+sign);
      oc1.setStatus(Particle::Status::Undecayed);

      //--- second outgoing central particle
      auto& oc2 = (*event_)[Particle::CentralSystem][1].get();
      oc2.setChargeSign(-sign);
      oc2.setStatus(Particle::Status::Undecayed);
    }

    //----- utilities

    double Process2to4::that() const {
      const double that1 = (q1() - pc(0)).mass2();
      const double that2 = (q2() - pc(1)).mass2();
      return 0.5 * (that1 + that2);
    }

    double Process2to4::uhat() const {
      const double uhat1 = (q1() - pc(1)).mass2();
      const double uhat2 = (q2() - pc(0)).mass2();
      return 0.5 * (uhat1 + uhat2);
    }

    ParametersDescription Process2to4::description() {
      auto desc = KTProcess::description();
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
