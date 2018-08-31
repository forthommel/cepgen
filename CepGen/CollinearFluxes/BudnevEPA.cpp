/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
 *                2003  Maarten Boonekamp, Tibor Kucs
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

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class BudnevEPA : public CollinearFlux {
  public:
    explicit BudnevEPA(const ParametersList& params) : CollinearFlux(params), q2range_(steer<Limits>("q2range")) {}

    bool fragmenting() const override final { return false; }
    pdgid_t partonPdgId() const override final { return PDG::photon; }

    static ParametersDescription description() {
      auto desc = CollinearFlux::description();
      desc.setDescription("Generic Budnev EPA flux");
      desc.add<Limits>("q2range", {0., 1.e4}).setDescription("kinematic range for the parton virtuality, in GeV^2");
      return desc;
    }

  protected:
    bool computeQ2min(double x, double& q2min) const {
      if (!x_range_.contains(x, true))
        return false;
      q2min = mass2() * x * x / (1. - x);
      return q2range_.contains(q2min);
    }

    const Limits q2range_;
  };

  class BudnevEPALepton final : public BudnevEPA {
  public:
    explicit BudnevEPALepton(const ParametersList& params)
        : BudnevEPA(params), ml2_(std::pow(PDG::get().mass(steer<int>("pdgId")), 2)) {
      CG_INFO("BudnevEPALepton") << "Budnev EPA for photon-from-lepton elastic limit (lepton: "
                                 << PDG::get().name(steer<int>("pdgId")) << ").\n\t "
                                 << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
    }

    static ParametersDescription description() {
      auto desc = BudnevEPA::description();
      desc.setDescription("Budnev EPA for lepton");
      desc.addAs<int, pdgid_t>("pdgId", PDG::electron).setDescription("lepton PDG id");
      return desc;
    }

    double fluxMX2(double x, double /*mf2*/) const override {
      double q2min;
      if (!computeQ2min(x, q2min))
        return 0.;
      return std::max(prefactor_ * (1. * 0.5 * ml2_ * x * (-1. / q2min + 1. / q2range_.max()) +
                                    (1. - x + 0.5 * x * x) / x * log(q2range_.max() / q2min)),
                      0.);
    }

  private:
    double mass2() const override { return ml2_; }
    const double ml2_;
  };

  class BudnevEPANucleon : public BudnevEPA {
  public:
    explicit BudnevEPANucleon(const ParametersList& params)
        : BudnevEPA(params),
          q2scale_(steer<double>("scale")),
          a_(steer<double>("a")),
          b_(steer<double>("b")),
          c_(steer<double>("c")) {}

    static ParametersDescription description() {
      auto desc = BudnevEPA::description();
      desc.add<double>("scale", 0.71);
      desc.add<double>("a", 7.16);
      desc.add<double>("b", -3.96);
      desc.add<double>("c", 0.028);
      return desc;
    }

    double fluxMX2(double x, double /*mf2*/) const override {
      double q2min;
      if (!computeQ2min(x, q2min))
        return 0.;
      return std::max(0., prefactor_ * (phi(x, q2range_.max() / q2scale_) - phi(x, q2min / q2scale_)) * (1 - x) / x);
    }

  protected:
    double phi(double x, double qq) const {
      if (x == 1.)
        return 0.;
      const double qq1 = 1 + qq, y = x * x / (1 - x);
      double f = (1 + a_ * y) * (-log(qq1 / qq) + 1 / qq1 + 1 / (2 * qq1 * qq1) + 1 / (3 * qq1 * qq1 * qq1));
      f += (1 - b_) * y / (4 * qq * qq1 * qq1 * qq1);
      f += c_ * (1 + y / 4) *
           (log((qq1 - b_) / qq1) + b_ / qq1 + b_ * b_ / (2 * qq1 * qq1) + b_ * b_ * b_ / (3 * qq1 * qq1 * qq1));
      return f;
    }

  private:
    const double q2scale_;
    const double a_, b_, c_;
  };

  struct BudnevEPAProton final : public BudnevEPANucleon {
    explicit BudnevEPAProton(const ParametersList& params) : BudnevEPANucleon(params) {
      CG_INFO("BudnevEPAProton") << "Budnev EPA for photon-from-proton elastic limit.\n\t"
                                 << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
    }

    double mass2() const override { return mp2_; }

    static ParametersDescription description() {
      auto desc = BudnevEPANucleon::description();
      desc.setDescription("Budnev EPA for proton");
      return desc;
    }
  };

  class BudnevEPAHI final : public BudnevEPANucleon {
  public:
    explicit BudnevEPAHI(const ParametersList& params)
        : BudnevEPANucleon(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          mass2_(hi_.mass() * hi_.mass()) {
      CG_INFO("BudnevEPAHI") << "Budnev EPA for photon-from-heavy ion elastic limit (HI: " << hi_ << ").\n\t"
                             << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      if (q2range_.max() < q2max_min_)
        throw CG_FATAL("BudnevEPAHI") << "Invalid maximal Q^2. Its value should not be below " << q2max_min_
                                      << " GeV^2. It is currently " << q2range_ << ".";
    }

    double mass2() const override { return mass2_; }

    static ParametersDescription description() {
      auto desc = BudnevEPANucleon::description();
      desc.setDescription("Budnev EPA for heavy ion");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb()).setDescription("type of heavy ion considered");
      desc.add<Limits>("q2range", {0., 1.e5}).setDescription("kinematic range for the parton virtuality, in GeV^2");
      return desc;
    }

    double fluxMX2(double x, double mf2) const override {
      const auto z = (unsigned short)hi_.Z;
      return z * BudnevEPANucleon::fluxMX2(x, mf2);
    }

  private:
    const HeavyIon hi_;
    const double mass2_;
    const double q2max_min_{1.e4};
  };
}  // namespace cepgen
REGISTER_FLUX("BudnevEPALepton", BudnevEPALepton);
REGISTER_FLUX("BudnevEPAHI", BudnevEPAHI);
REGISTER_FLUX("BudnevEPAProton", BudnevEPAProton);
