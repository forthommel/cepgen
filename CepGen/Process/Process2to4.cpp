/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  namespace proc {
    Process2to4::Process2to4(const ParametersList& params, pdgid_t cs_id)
        : FactorisedProcess(params, {cs_id, cs_id}), cs_prop_(PDG::get()(cs_id)), single_limits_(params) {}

    void Process2to4::setCuts(const cuts::Central& single) { single_limits_ = single; }

    void Process2to4::prepareFactorisedPhaseSpace() {
      if (cs_prop_.pdgid == PDG::invalid)  // ensure the central particles properties are correctly initialised
        cs_prop_ = PDG::get()(steer<ParticleProperties>("pair").pdgid);

      const auto lim_pt_diff = kinematics().cuts().central.pt_diff.truncate(Limits{0., 500.});
      defineVariable(
          m_pt_diff_, Mapping::linear, lim_pt_diff, "pt_diff", "Final state particles transverse momentum difference");

      const auto lim_phi_diff = kinematics().cuts().central.phi_diff.truncate(Limits{0., 2. * M_PI});
      defineVariable(m_phi_pt_diff_,
                     Mapping::linear,
                     lim_phi_diff,
                     "phi_pt_diff",
                     "Final state particles azimuthal angle difference");

      prepareProcessKinematics();

      lim_rap_ = kinematics().cuts().central.rapidity_single.truncate(Limits{-6., 6.});
      inv_sqrts_ = inverseSqrtS();
    }

    double Process2to4::computeFactorisedMatrixElement() {
      double jacob = 1.;
      {
        const auto xprod = x1() * x2();
        if (!utils::positive(xprod))
          return 0.;
        const auto qt_sum = (q1() + q2()).transverse();  // two-parton system
        const auto pt_diff = Momentum::fromPtEtaPhiE(m_pt_diff_, 0., m_phi_pt_diff_);
        const auto pt_c1 = 0.5 * (qt_sum + pt_diff), pt_c2 = 0.5 * (qt_sum - pt_diff);
        const auto p1t = pt_c1.pt(), p2t = pt_c2.pt();
        // apply user cuts on central system
        if (!kinematics().cuts().central.pt_single.contains(p1t) || !single_limits_.pt_single.contains(p1t))
          return 0.;
        if (!kinematics().cuts().central.pt_single.contains(p2t) || !single_limits_.pt_single.contains(p2t))
          return 0.;
        if (!kinematics().cuts().central.pt_diff.contains(std::fabs(p1t - p2t)))  // transverse momentum difference
          return 0.;
        // compute central particles rapidity
        const auto amt1 = utils::fastHypot(cs_prop_.mass, p1t) * inv_sqrts_,
                   amt2 = utils::fastHypot(cs_prop_.mass, p2t) * inv_sqrts_, amt2_diff = amt1 * amt1 - amt2 * amt2;
        if (!utils::positive(amt1 * amt2))
          return 0.;
        const auto xfrac = std::sqrt(std::pow(amt2_diff + xprod, 2) - 4. * amt1 * amt1 * xprod) + xprod;
        const auto y_c1 = +std::log(0.5 * (xfrac + amt2_diff) / amt1 / x2()),
                   y_c2 = -std::log(0.5 * (xfrac - amt2_diff) / amt2 / x1());
        if (!lim_rap_.contains(y_c1) || !lim_rap_.contains(y_c2))  // single particle rapidity
          return 0.;
        const auto y_diff = std::fabs(y_c1 - y_c2);
        if (!kinematics().cuts().central.rapidity_diff.contains(y_diff))  // rapidity distance
          return 0.;
        jacob = amt1 * amt2 * std::sinh(y_diff);
        if (!utils::positive(jacob))
          return 0.;
        // compute the four-momenta of the outgoing central particles
        pc(0) = Momentum::fromPtYPhiM(p1t, y_c1, pt_c1.phi(), cs_prop_.mass);
        pc(1) = Momentum::fromPtYPhiM(p2t, y_c2, pt_c2.phi(), cs_prop_.mass);

        CG_DEBUG_LOOP("Process2to4:computeFactorisedMatrixElement").log([&](auto& log) {
          const auto x1_comp = amt1 * exp(+y_c1) + amt2 * exp(+y_c2);
          const auto x2_comp = amt1 * exp(-y_c1) + amt2 * exp(-y_c2);
          log << "Closure test for x_i computation:\n"
              << "\tx1(gen)=" << x1() << ", x1(comp)=" << x1_comp << ", diff=" << (x1_comp - x1()) << "\n"
              << "\tx2(gen)=" << x2() << ", x2(comp)=" << x2_comp << ", diff=" << (x2_comp - x2()) << ".";
        });
      }

      //--- window in central system invariant mass
      const auto invm = (pc(0) + pc(1)).mass();
      if (!kinematics().cuts().central.mass_sum.contains(invm))
        return 0.;

      //--- additional conditions for energy-momentum conservation
      if (!kinematics().incomingBeams().positive().elastic() && std::sqrt(x2() * s() - invm - q2().p2()) <= mX())
        return 0.;
      if (!kinematics().incomingBeams().negative().elastic() && std::sqrt(x1() * s() - invm - q1().p2()) <= mY())
        return 0.;

      //--- four-momenta of the outgoing protons (or remnants)

      const auto px_p = (1. - x1()) * pA().p() * M_SQRT2, px_m = (mX2() + q1().p2()) * 0.5 / px_p;
      const auto py_m = (1. - x2()) * pB().p() * M_SQRT2, py_p = (mY2() + q2().p2()) * 0.5 / py_m;
      CG_DEBUG_LOOP("2to4:pxy") << "px+ = " << px_p << " / px- = " << px_m << "\n\t"
                                << "py+ = " << py_p << " / py- = " << py_m << ".";

      pX() = -Momentum(q1()).setPz((px_p - px_m) * M_SQRT1_2).setEnergy((px_p + px_m) * M_SQRT1_2);
      pY() = -Momentum(q2()).setPz((py_p - py_m) * M_SQRT1_2).setEnergy((py_p + py_m) * M_SQRT1_2);

      CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << pX() << ", mass = " << pX().mass() << "\n\t"
                                     << "Second remnant: " << pY() << ", mass = " << pY().mass() << ".";

      if (std::fabs(pX().mass2() - mX2()) > NUM_LIMITS) {
        CG_WARNING("2to4:px") << "Invalid X system squared mass: " << pX().mass2() << "/" << mX2() << ".";
        return 0.;
      }
      if (std::fabs(pY().mass2() - mY2()) > NUM_LIMITS) {
        CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << pY().mass2() << "/" << mY2() << ".";
        return 0.;
      }

      //--- four-momenta of the intermediate partons
      const double norm = 1. / wCM() / wCM() / s(), prefac = 0.5 * wCM() * sqrtS();
      {  // positive-z incoming parton collinear kinematics
        const double tau1 = norm * q1().p2() / x1() / x1();
        q1().setPz(+prefac * x1() * (1. - tau1)).setEnergy(+prefac * x1() * (1. + tau1));
      }
      {  // negative-z incoming parton collinear kinematics
        const double tau2 = norm * q2().p2() / x2() / x2();
        q2().setPz(-prefac * x2() * (1. - tau2)).setEnergy(+prefac * x2() * (1. + tau2));
      }

      CG_DEBUG_LOOP("2to4:partons") << "Squared c.m. energy = " << s() << " GeV^2\n\t"
                                    << "First parton: " << q1() << ", mass2 = " << q1().mass2() << ", x1 = " << x1()
                                    << ", p = " << q1().p() << "\n\t"
                                    << "Second parton: " << q2() << ", mass2 = " << q2().mass2() << ", x2 = " << x2()
                                    << ", p = " << q2().p() << ".";

      if (const auto amat2 = computeCentralMatrixElement(); utils::positive(amat2))
        return amat2 * prefactor_ * m_pt_diff_ / jacob;
      return 0.;  // skip computing the prefactors if invalid
    }

    void Process2to4::fillCentralParticlesKinematics() {
      const short sign = rnd_gen_->uniformInt(0, 1) == 1 ? 1 : -1;  // randomise the charge of outgoing system
      event()[Particle::CentralSystem][0].get().setChargeSign(+sign).setStatus(Particle::Status::Undecayed);
      event()[Particle::CentralSystem][1].get().setChargeSign(-sign).setStatus(Particle::Status::Undecayed);
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
