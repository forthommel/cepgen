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

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PhaseSpaceGenerator.h"
#include "CepGen/Process/KTProcess.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  void PhaseSpaceGenerator2to3::initialise() {}

  bool PhaseSpaceGenerator2to3::generate() { return true; }

  void PhaseSpaceGenerator2to4::initialise() {
    proc_.defineVariable(y_c1_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.rapidity_single,
                         {-6., 6.},
                         "First outgoing particle rapidity");
    proc_.defineVariable(y_c2_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.rapidity_single,
                         {-6., 6.},
                         "Second outgoing particle rapidity");
    proc_.defineVariable(pt_diff_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.pt_diff,
                         {0., 500.},
                         "Final state particles transverse momentum difference");
    proc_.defineVariable(phi_pt_diff_,
                         proc::Process::Mapping::linear,
                         proc_.kinematics().cuts().central.phi_diff,
                         {0., 2. * M_PI},
                         "Final state particles azimuthal angle difference");
  }

  bool PhaseSpaceGenerator2to4::generate() {
    //--- transverse kinematics of outgoing central system
    const auto pt_diff = Momentum::fromPtEtaPhiE(pt_diff_, 0., phi_pt_diff_);
    if (fabs(pt_diff.pt() - pt_diff_) > proc_.NUM_LIMITS)
      throw CG_FATAL("Process2to4") << "|dpt|=" << pt_diff_ << " != dpt.pt()=" << pt_diff.pt() << ", dpt=" << pt_diff
                                    << ".";

    //--- two-parton system (in transverse plane)
    const auto qt_sum = proc_.q1() + proc_.q2();

    const auto pt_c1 = 0.5 * (qt_sum + pt_diff);
    const auto pt_c2 = 0.5 * (qt_sum - pt_diff);
    const double p1t = pt_c1.pt(), p2t = pt_c2.pt();

    CG_DEBUG_LOOP("2to4:me") << "diff(pt) = " << pt_diff << "\n\t"
                             << "p(1/2)x = " << pt_c1.px() << " / " << pt_c2.px() << "\n\t"
                             << "p(1/2)y = " << pt_c1.py() << " / " << pt_c2.py() << "\n\t"
                             << "p(1/2)t = " << p1t << " / " << p2t;

    //--- window in rapidity distance
    if (!proc_.kinematics().cuts().central.rapidity_diff.contains(fabs(y_c1_ - y_c2_)))
      return false;

    //--- apply the pt cut already at this stage (remains unchanged)
    if (!proc_.kinematics().cuts().central.pt_single.contains(p1t))
      return false;
    if (!proc_.kinematics().cuts().central.pt_single.contains(p2t))
      return false;
    if (!single_limits_.pt_single.contains(p1t))
      return false;
    if (!single_limits_.pt_single.contains(p2t))
      return false;

    //--- window in transverse momentum difference
    if (!proc_.kinematics().cuts().central.pt_diff.contains(fabs(p1t - p2t)))
      return false;

    //--- transverse mass for the two central particles
    amt1_ = std::hypot(p1t, cs_prop_.mass);
    amt2_ = std::hypot(p2t, cs_prop_.mass);

    //--- window in central system invariant mass
    const double invm = sqrt(amt1_ * amt1_ + amt2_ * amt2_ + 2. * amt1_ * amt2_ * cosh(y_c1_ - y_c2_) - qt_sum.pt2());
    if (!proc_.kinematics().cuts().central.mass_sum.contains(invm))
      return false;

    //--- auxiliary quantities

    const double alpha1 = amt1_ / sqs_ * exp(y_c1_), beta1 = amt1_ / sqs_ * exp(-y_c1_);
    const double alpha2 = amt2_ / sqs_ * exp(y_c2_), beta2 = amt2_ / sqs_ * exp(-y_c2_);
    x1_ = alpha1 + alpha2;
    x2_ = beta1 + beta2;

    CG_DEBUG_LOOP("2to4:sudakov") << "Sudakov parameters:\n\t"
                                  << "  alpha(1/2) = " << alpha1 << " / " << alpha2 << "\n\t"
                                  << "   beta(1/2) = " << beta1 << " / " << beta2 << ".";

    //--- sanity check for x_i values
    if (!x_limits_.contains(x1_) || !x_limits_.contains(x2_))
      return false;

    //--- additional conditions for energy-momentum conservation

    const double s1_eff = x1_ * s_ - qt1_ * qt1_, s2_eff = x2_ * s_ - qt2_ * qt2_;

    CG_DEBUG_LOOP("2to4:central") << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
                                  << "central system invariant mass = " << invm << " GeV";

    if (proc_.kinematics().incomingBeams().positive().fragmented() && (sqrt(s2_eff) <= sqrt(mX2_) + invm))
      return false;
    if (proc_.kinematics().incomingBeams().negative().fragmented() && (sqrt(s1_eff) <= sqrt(mY2_) + invm))
      return false;

    //--- four-momenta of the outgoing protons (or remnants)

    const double px_plus = (1. - x1_) * proc_.pA().p() * M_SQRT2;
    const double py_minus = (1. - x2_) * proc_.pB().p() * M_SQRT2;
    const double px_minus = (mX2_ + qt1_ * qt1_) * 0.5 / px_plus;
    const double py_plus = (mY2_ + qt2_ * qt2_) * 0.5 / py_minus;
    // warning! sign of pz??

    CG_DEBUG_LOOP("2to4:pxy") << "px± = " << px_plus << " / " << px_minus << "\n\t"
                              << "py± = " << py_plus << " / " << py_minus << ".";

    proc_.pX() =
        -Momentum(proc_.q1()).setPz((px_plus - px_minus) * M_SQRT1_2).setEnergy((px_plus + px_minus) * M_SQRT1_2);
    proc_.pY() =
        -Momentum(proc_.q2()).setPz((py_plus - py_minus) * M_SQRT1_2).setEnergy((py_plus + py_minus) * M_SQRT1_2);

    CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << proc_.pX() << ", mass = " << proc_.pX().mass() << "\n\t"
                                   << "Second remnant: " << proc_.pY() << ", mass = " << proc_.pY().mass() << ".";

    if (fabs(proc_.pX().mass2() - mX2_) > proc_.NUM_LIMITS) {
      CG_WARNING("2to4:px") << "Invalid X system squared mass: " << proc_.pX().mass2() << "/" << mX2_ << ".";
      return false;
    }
    if (fabs(proc_.pY().mass2() - mY2_) > proc_.NUM_LIMITS) {
      CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << proc_.pY().mass2() << "/" << mY2_ << ".";
      return false;
    }

    //--- four-momenta of the intermediate partons
    const double norm = 1. / ww_ / ww_ / s_;
    const double tau1 = norm * qt1_ * qt1_ / x1_ / x1_;
    proc_.q1().setPz(+0.5 * x1_ * ww_ * sqs_ * (1. - tau1)).setEnergy(+0.5 * x1_ * ww_ * sqs_ * (1. + tau1));

    const double tau2 = norm * qt2_ * qt2_ / x2_ / x2_;
    proc_.q2().setPz(-0.5 * x2_ * ww_ * sqs_ * (1. - tau2)).setEnergy(+0.5 * x2_ * ww_ * sqs_ * (1. + tau2));

    CG_DEBUG_LOOP("2to4:partons") << "First parton:  " << proc_.q1() << ", mass2 = " << proc_.q1().mass2() << "\n\t"
                                  << "Second parton: " << proc_.q2() << ", mass2 = " << proc_.q2().mass2() << ".";

    //--- four-momenta of the outgoing central particles

    proc_.pc(0) = (pt_c1 + alpha1 * proc_.pA() + beta1 * proc_.pB())
                      .setEnergy(alpha1 * proc_.pA().energy() + beta1 * proc_.pB().energy());
    proc_.pc(1) = (pt_c2 + alpha2 * proc_.pA() + beta2 * proc_.pB())
                      .setEnergy(alpha2 * proc_.pA().energy() + beta2 * proc_.pB().energy());

    CG_DEBUG_LOOP("2to4:central") << "First central particle:  " << proc_.pc(0) << ", mass = " << proc_.pc(0).mass()
                                  << "\n\t"
                                  << "Second central particle: " << proc_.pc(1) << ", mass = " << proc_.pc(1).mass()
                                  << ".";
    return true;
  }
}  // namespace cepgen
