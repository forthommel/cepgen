/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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
#include "CepGen/Process/Process2to4.h"

namespace cepgen {
  namespace proc {
    const Limits Process2to4::x_limits_{0., 1.};

    Process2to4::Process2to4(const ParametersList& params, pdgid_t cs_id)
        : KTProcess(params, {cs_id, cs_id}), cs_prop_(PDG::get()(cs_id)), single_limits_(params), rnd_sign_(0, 1) {}

    void Process2to4::setCuts(const cuts::Central& single) { single_limits_ = single; }

    void Process2to4::preparePhaseSpace() {
      if (cs_prop_.pdgid == PDG::invalid)  // ensure the central particles properties are correctly initialised
        cs_prop_ = PDG::get()(steer<ParticleProperties>("pair").pdgid);

      const auto lim_rap = kinematics().cuts().central.rapidity_single.truncate(Limits{-6., 6.});
      defineVariable(m_y_c1_, Mapping::linear, lim_rap, "First outgoing particle rapidity");
      defineVariable(m_y_c2_, Mapping::linear, lim_rap, "Second outgoing particle rapidity");

      const auto lim_pt_diff = kinematics().cuts().central.pt_diff.truncate(Limits{0., 500.});
      defineVariable(m_pt_diff_, Mapping::linear, lim_pt_diff, "Final state particles transverse momentum difference");

      const auto lim_phi_diff = kinematics().cuts().central.phi_diff.truncate(Limits{0., 2. * M_PI});
      defineVariable(m_phi_pt_diff_, Mapping::linear, lim_phi_diff, "Final state particles azimuthal angle difference");

      ww_ = 0.5 * (1. + std::sqrt(1. - 4. * mA() * mB() / s()));

      prepareProcessKinematics();
    }

    double Process2to4::computeKTFactorisedMatrixElement() {
      //--- transverse kinematics of initial partons
      const auto qt_1 = Momentum::fromPtEtaPhiE(m_qt1_, 0., m_phi_qt1_),
                 qt_2 = Momentum::fromPtEtaPhiE(m_qt2_, 0., m_phi_qt2_), qt_sum = qt_1 + qt_2;

      //--- transverse kinematics of outgoing central system
      const auto pt_diff = Momentum::fromPtEtaPhiE(m_pt_diff_, 0., m_phi_pt_diff_);

      {  //--- define the central system from the computed quantities
        const auto pt_c1 = 0.5 * (qt_sum + pt_diff);
        const auto pt_c2 = 0.5 * (qt_sum - pt_diff);
        const double p1t = pt_c1.pt(), p2t = pt_c2.pt();
        //--- cuts on central system
        if (!kinematics().cuts().central.rapidity_diff.contains(fabs(m_y_c1_ - m_y_c2_)))  // rapidity distance
          return 0.;
        if (!kinematics().cuts().central.pt_single.contains(p1t) || !single_limits_.pt_single.contains(p1t))
          return 0.;
        if (!kinematics().cuts().central.pt_single.contains(p2t) || !single_limits_.pt_single.contains(p2t))
          return 0.;
        if (!kinematics().cuts().central.pt_diff.contains(fabs(p1t - p2t)))  // transverse momentum difference
          return 0.;
        //--- four-momenta of the outgoing central particles
        pc(0) = Momentum::fromPtYPhiM(p1t, m_y_c1_, pt_c1.phi(), cs_prop_.mass);
        pc(1) = Momentum::fromPtYPhiM(p2t, m_y_c2_, pt_c2.phi(), cs_prop_.mass);
      }

      //--- window in central system invariant mass
      const auto invm = (pc(0) + pc(1)).mass();
      if (!kinematics().cuts().central.mass_sum.contains(invm))
        return 0.;

      //--- auxiliary quantities

      //--- transverse mass for the two central particles
      const auto amt1 = pc(0).massT(), amt2 = pc(1).massT();

      x1_ = amt1 / sqrtS() * exp(+m_y_c1_) + amt2 / sqrtS() * exp(+m_y_c2_);
      x2_ = amt1 / sqrtS() * exp(-m_y_c1_) + amt2 / sqrtS() * exp(-m_y_c2_);

      //--- sanity check for x_i values
      if (!x_limits_.contains(x1_) || !x_limits_.contains(x2_))
        return 0.;

      //--- additional conditions for energy-momentum conservation

      if (!kinematics().incomingBeams().positive().elastic() && std::sqrt(x2_ * s() - invm - m_qt2_ * m_qt2_) <= mX())
        return 0.;
      if (!kinematics().incomingBeams().negative().elastic() && std::sqrt(x1_ * s() - invm - m_qt1_ * m_qt1_) <= mY())
        return 0.;

      //--- four-momenta of the outgoing protons (or remnants)

      const double px_plus = (1. - x1_) * pA().p() * M_SQRT2;
      const double py_minus = (1. - x2_) * pB().p() * M_SQRT2;
      const double px_minus = (mX2() + m_qt1_ * m_qt1_) * 0.5 / px_plus;
      const double py_plus = (mY2() + m_qt2_ * m_qt2_) * 0.5 / py_minus;
      // warning! sign of pz??

      CG_DEBUG_LOOP("2to4:pxy") << "px± = " << px_plus << " / " << px_minus << "\n\t"
                                << "py± = " << py_plus << " / " << py_minus << ".";

      pX() = -Momentum(qt_1).setPz((px_plus - px_minus) * M_SQRT1_2).setEnergy((px_plus + px_minus) * M_SQRT1_2);
      pY() = -Momentum(qt_2).setPz((py_plus - py_minus) * M_SQRT1_2).setEnergy((py_plus + py_minus) * M_SQRT1_2);

      CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << pX() << ", mass = " << pX().mass() << "\n\t"
                                     << "Second remnant: " << pY() << ", mass = " << pY().mass() << ".";

      if (fabs(pX().mass2() - mX2()) > NUM_LIMITS) {
        CG_WARNING("2to4:px") << "Invalid X system squared mass: " << pX().mass2() << "/" << mX2() << ".";
        return 0.;
      }
      if (fabs(pY().mass2() - mY2()) > NUM_LIMITS) {
        CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << pY().mass2() << "/" << mY2() << ".";
        return 0.;
      }

      //--- four-momenta of the intermediate partons

      const double norm = 1. / ww_ / ww_ / s();
      const double tau1 = norm * m_qt1_ * m_qt1_ / x1_ / x1_;
      q1() = Momentum(qt_1)
                 .setPz(+0.5 * x1_ * ww_ * sqrtS() * (1. - tau1))
                 .setEnergy(+0.5 * x1_ * ww_ * sqrtS() * (1. + tau1));

      const double tau2 = norm * m_qt2_ * m_qt2_ / x2_ / x2_;
      q2() = Momentum(qt_2)
                 .setPz(-0.5 * x2_ * ww_ * sqrtS() * (1. - tau2))
                 .setEnergy(+0.5 * x2_ * ww_ * sqrtS() * (1. + tau2));

      CG_DEBUG_LOOP("2to4:partons") << "First parton:  " << q1() << ", mass2 = " << q1().mass2() << "\n\t"
                                    << "Second parton: " << q2() << ", mass2 = " << q2().mass2() << ".";

      //--- compute the central 2-to-2 matrix element

      const double amat2 = computeCentralMatrixElement();
      if (amat2 <= 0.)  // skip computing the fluxes if no contribution
        return 0.;

      // factor 1/4 from jacobian of transformations
      return amat2 * std::pow(4. * M_PI * x1_ * x2_ * s(), -2) * 0.25 * constants::GEVM2_TO_PB *
             (m_pt_diff_ * m_qt1_ * m_qt2_);
    }

    void Process2to4::fillCentralParticlesKinematics() {
      //--- randomise the charge of outgoing system
      const short sign = rnd_sign_(rnd_gen_) == 1 ? 1 : -1;

      //--- first outgoing central particle
      auto& oc1 = event()[Particle::CentralSystem][0].get();
      oc1.setChargeSign(+sign);
      oc1.setStatus(Particle::Status::Undecayed);

      //--- second outgoing central particle
      auto& oc2 = event()[Particle::CentralSystem][1].get();
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
  }  // namespace proc
}  // namespace cepgen
