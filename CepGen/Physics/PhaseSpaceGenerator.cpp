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
  const Limits PhaseSpaceGenerator::x_limits_{0., 1.};

  PhaseSpaceGenerator::PhaseSpaceGenerator(proc::KTProcess& proc)
      : proc_(proc), ww_(0.5 * (1. + sqrt(1. - 4. * sqrt(proc_.mA2_ * proc_.mB2_) / proc_.s_))) {}

  void PhaseSpaceGenerator2to3::initialise() {
    proc_.defineVariable(y_c_,
                         proc::Process::Mapping::linear,
                         proc_.kin_.cuts().central.rapidity_single,
                         {-6., 6.},
                         "Ooutgoing particle rapidity");
    proc_.defineVariable(pt_c_,
                         proc::Process::Mapping::linear,
                         proc_.kin_.cuts().central.pt_single,
                         {0., 500.},
                         "Final state particle transverse momentum");
  }

  bool PhaseSpaceGenerator2to3::generate() {
    //--- transverse kinematics of initial partons
    const auto qt_1 = Momentum::fromPtEtaPhiE(proc_.qt1_, 0., proc_.phi_qt1_);
    if (fabs(qt_1.pt() - proc_.qt1_) > proc_.NUM_LIMITS)
      throw CG_FATAL("Process2to3") << "|qt1|=" << proc_.qt1_ << " != qt1.pt()=" << qt_1.pt() << ", qt1=" << qt_1
                                    << ".";

    const auto qt_2 = Momentum::fromPtEtaPhiE(proc_.qt2_, 0., proc_.phi_qt2_);
    if (fabs(qt_2.pt() - proc_.qt2_) > proc_.NUM_LIMITS)
      throw CG_FATAL("Process2to3") << "|qt2|=" << proc_.qt2_ << " != qt2.pt()=" << qt_2.pt() << ", qt2=" << qt_2
                                    << ".";

    //--- two-parton system (in transverse plane)
    const auto qt_sum = qt_1 + qt_2;

    CG_DEBUG_LOOP("Process2to3:me") << "q(1/2)x = " << qt_1.px() << " / " << qt_2.px() << "\n\t"
                                    << "q(1/2)y = " << qt_1.py() << " / " << qt_2.py() << "\n\t"
                                    << "sum(qt) = " << qt_sum;

    //--- transverse kinematics of outgoing central system
    const auto pt_c = Momentum::fromPtEtaPhiE(pt_c_, 0., qt_sum.phi());
    const double pt = pt_c.pt();

    CG_DEBUG_LOOP("Process2to3:me") << "pc = " << pt_c << ", ptc = " << pt << ".";

    // window in rapidity
    if (!proc_.kin_.cuts().central.rapidity_single.contains(y_c_))
      return 0.;

    // apply the pt cut already at this stage (remains unchanged)
    if (!proc_.kin_.cuts().central.pt_single.contains(pt))
      return 0.;
    /*if (!single_limits_.pt_single.contains(pt))
      return 0.;*/ //FIXME

    // transverse mass for the central particle
    proc_.mt_.at(0) = std::hypot(pt, PDG::get()(proc_.produced_parts_.at(0)).mass);

    // window in central system invariant mass
    const double invm = proc_.mt_.at(0) * sqrt(0.5 + 0.5 * cosh(y_c_) - qt_sum.pt2());
    if (!proc_.kin_.cuts().central.mass_sum.contains(invm))
      return 0.;

    // auxiliary quantities
    const double x1 = proc_.mt_.at(0) / proc_.sqs_ * exp(+y_c_), x2 = proc_.mt_.at(0) / proc_.sqs_ * exp(-y_c_);

    CG_DEBUG_LOOP("Process2to3:sudakov") << "Sudakov parameters:\n\t"
                                         << "  x(1/2) = " << std::make_pair(x1, x2) << ".";

    //--- sanity check for x_i values
    if (!x_limits_.contains(x1) || !x_limits_.contains(x2))
      return 0.;

    //--- additional conditions for energy-momentum conservation

    const double s1_eff = x1 * proc_.s_ - proc_.qt1_ * proc_.qt1_, s2_eff = x2 * proc_.s_ - proc_.qt2_ * proc_.qt2_;

    CG_DEBUG_LOOP("Process2to3:central") << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
                                         << "central system invariant mass = " << invm << " GeV";

    if (proc_.kin_.incomingBeams().positive().fragmented() && (sqrt(s2_eff) <= sqrt(proc_.mX2_) + invm))
      return 0.;
    if (proc_.kin_.incomingBeams().negative().fragmented() && (sqrt(s1_eff) <= sqrt(proc_.mY2_) + invm))
      return 0.;

    // four-momenta of the outgoing protons (or remnants)
    const double px_plus = (1. - x1) * proc_.pA().p() * M_SQRT2;
    const double py_minus = (1. - x2) * proc_.pB().p() * M_SQRT2;
    const double px_minus = (proc_.mX2_ + proc_.qt1_ * proc_.qt1_) * 0.5 / px_plus;
    const double py_plus = (proc_.mY2_ + proc_.qt2_ * proc_.qt2_) * 0.5 / py_minus;
    proc_.pX() = -Momentum(qt_1).setPz((px_plus - px_minus) * M_SQRT1_2).setEnergy((px_plus + px_minus) * M_SQRT1_2);
    proc_.pY() = -Momentum(qt_2).setPz((py_plus - py_minus) * M_SQRT1_2).setEnergy((py_plus + py_minus) * M_SQRT1_2);

    if (fabs(proc_.pX().mass2() - proc_.mX2_) > proc_.NUM_LIMITS) {
      CG_WARNING("Process2to3:px") << "Invalid X system squared mass: " << proc_.pX().mass2() << "/" << proc_.mX2_
                                   << ".";
      return 0.;
    }
    if (fabs(proc_.pY().mass2() - proc_.mY2_) > proc_.NUM_LIMITS) {
      CG_WARNING("Process2to3:py") << "Invalid Y system squared mass: " << proc_.pY().mass2() << "/" << proc_.mY2_
                                   << ".";
      return 0.;
    }

    // four-momenta of the intermediate partons
    const double norm = 1. / ww_ / ww_ / proc_.s_;
    const double tau1 = norm * proc_.qt1_ * proc_.qt1_ / x1 / x1;
    proc_.x1_ = x1;
    proc_.q1() = Momentum(qt_1)
                     .setPz(+0.5 * x1 * ww_ * proc_.sqs_ * (1. - tau1))
                     .setEnergy(+0.5 * x1 * ww_ * proc_.sqs_ * (1. + tau1));

    const double tau2 = norm * proc_.qt2_ * proc_.qt2_ / x2 / x2;
    proc_.x2_ = x2;
    proc_.q2() = Momentum(qt_2)
                     .setPz(-0.5 * x2 * ww_ * proc_.sqs_ * (1. - tau2))
                     .setEnergy(+0.5 * x2 * ww_ * proc_.sqs_ * (1. + tau2));

    // four-momenta of the outgoing central particles
    proc_.pc(0) = (pt_c + proc_.x1_ * proc_.pA() + proc_.x2_ * proc_.pB())
                      .setEnergy(proc_.x1_ * proc_.pA().energy() + proc_.x2_ * proc_.pB().energy());
    return true;
  }

  void PhaseSpaceGenerator2to4::initialise() {
    proc_.defineVariable(y_c1_,
                         proc::Process::Mapping::linear,
                         proc_.kin_.cuts().central.rapidity_single,
                         {-6., 6.},
                         "First outgoing particle rapidity");
    proc_.defineVariable(y_c2_,
                         proc::Process::Mapping::linear,
                         proc_.kin_.cuts().central.rapidity_single,
                         {-6., 6.},
                         "Second outgoing particle rapidity");
    proc_.defineVariable(pt_diff_,
                         proc::Process::Mapping::linear,
                         proc_.kin_.cuts().central.pt_diff,
                         {0., 500.},
                         "Final state particles transverse momentum difference");
    proc_.defineVariable(phi_pt_diff_,
                         proc::Process::Mapping::linear,
                         proc_.kin_.cuts().central.phi_diff,
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
    if (!proc_.kin_.cuts().central.rapidity_diff.contains(fabs(y_c1_ - y_c2_)))
      return false;

    //--- apply the pt cut already at this stage (remains unchanged)
    if (!proc_.kin_.cuts().central.pt_single.contains(p1t))
      return false;
    if (!proc_.kin_.cuts().central.pt_single.contains(p2t))
      return false;
    /*if (!single_limits_.pt_single.contains(p1t))
      return false;
    if (!single_limits_.pt_single.contains(p2t))
      return false;*/ //FIXME

    //--- window in transverse momentum difference
    if (!proc_.kin_.cuts().central.pt_diff.contains(fabs(p1t - p2t)))
      return false;

    //--- transverse mass for the two central particles
    proc_.mt_.at(0) = std::hypot(p1t, PDG::get()(proc_.produced_parts_.at(0)).mass);
    proc_.mt_.at(1) = std::hypot(p2t, PDG::get()(proc_.produced_parts_.at(1)).mass);

    //--- window in central system invariant mass
    const double invm = sqrt(proc_.mt_.at(0) * proc_.mt_.at(0) + proc_.mt_.at(1) * proc_.mt_.at(1) +
                             2. * proc_.mt_.at(0) * proc_.mt_.at(1) * cosh(y_c1_ - y_c2_) - qt_sum.pt2());
    if (!proc_.kin_.cuts().central.mass_sum.contains(invm))
      return false;

    //--- auxiliary quantities

    const double alpha1 = proc_.mt_.at(0) / proc_.sqs_ * exp(y_c1_), beta1 = proc_.mt_.at(0) / proc_.sqs_ * exp(-y_c1_);
    const double alpha2 = proc_.mt_.at(1) / proc_.sqs_ * exp(y_c2_), beta2 = proc_.mt_.at(1) / proc_.sqs_ * exp(-y_c2_);
    const double x1 = alpha1 + alpha2, x2 = beta1 + beta2;

    CG_DEBUG_LOOP("2to4:sudakov") << "Sudakov parameters:\n\t"
                                  << "  alpha(1/2) = " << alpha1 << " / " << alpha2 << "\n\t"
                                  << "   beta(1/2) = " << beta1 << " / " << beta2 << ".";

    //--- sanity check for x_i values
    if (!x_limits_.contains(x1) || !x_limits_.contains(x2))
      return false;

    //--- additional conditions for energy-momentum conservation

    const double s1_eff = x1 * proc_.s_ - proc_.qt1_ * proc_.qt1_, s2_eff = x2 * proc_.s_ - proc_.qt2_ * proc_.qt2_;

    CG_DEBUG_LOOP("2to4:central") << "s(1/2)_eff = " << s1_eff << " / " << s2_eff << " GeV^2\n\t"
                                  << "central system invariant mass = " << invm << " GeV";

    if (proc_.kin_.incomingBeams().positive().fragmented() && (sqrt(s2_eff) <= sqrt(proc_.mX2_) + invm))
      return false;
    if (proc_.kin_.incomingBeams().negative().fragmented() && (sqrt(s1_eff) <= sqrt(proc_.mY2_) + invm))
      return false;

    //--- four-momenta of the outgoing protons (or remnants)

    const double px_plus = (1. - x1) * proc_.pA().p() * M_SQRT2;
    const double py_minus = (1. - x2) * proc_.pB().p() * M_SQRT2;
    const double px_minus = (proc_.mX2_ + proc_.qt1_ * proc_.qt1_) * 0.5 / px_plus;
    const double py_plus = (proc_.mY2_ + proc_.qt2_ * proc_.qt2_) * 0.5 / py_minus;
    // warning! sign of pz??

    CG_DEBUG_LOOP("2to4:pxy") << "px± = " << px_plus << " / " << px_minus << "\n\t"
                              << "py± = " << py_plus << " / " << py_minus << ".";

    proc_.pX() =
        -Momentum(proc_.q1()).setPz((px_plus - px_minus) * M_SQRT1_2).setEnergy((px_plus + px_minus) * M_SQRT1_2);
    proc_.pY() =
        -Momentum(proc_.q2()).setPz((py_plus - py_minus) * M_SQRT1_2).setEnergy((py_plus + py_minus) * M_SQRT1_2);

    CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << proc_.pX() << ", mass = " << proc_.pX().mass() << "\n\t"
                                   << "Second remnant: " << proc_.pY() << ", mass = " << proc_.pY().mass() << ".";

    if (fabs(proc_.pX().mass2() - proc_.mX2_) > proc_.NUM_LIMITS) {
      CG_WARNING("2to4:px") << "Invalid X system squared mass: " << proc_.pX().mass2() << "/" << proc_.mX2_ << ".";
      return false;
    }
    if (fabs(proc_.pY().mass2() - proc_.mY2_) > proc_.NUM_LIMITS) {
      CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << proc_.pY().mass2() << "/" << proc_.mY2_ << ".";
      return false;
    }

    //--- four-momenta of the intermediate partons
    const double norm = 1. / ww_ / ww_ / proc_.s_;
    const double tau1 = norm * proc_.qt1_ * proc_.qt1_ / x1 / x1;
    proc_.x1_ = x1;
    proc_.q1().setPz(+0.5 * x1 * ww_ * proc_.sqs_ * (1. - tau1)).setEnergy(+0.5 * x1 * ww_ * proc_.sqs_ * (1. + tau1));

    const double tau2 = norm * proc_.qt2_ * proc_.qt2_ / proc_.x2_ / proc_.x2_;
    proc_.x2_ = x2;
    proc_.q2().setPz(-0.5 * x2 * ww_ * proc_.sqs_ * (1. - tau2)).setEnergy(+0.5 * x2 * ww_ * proc_.sqs_ * (1. + tau2));

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
