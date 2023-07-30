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

#include "CepGen/CollinearFluxes/IntegratedPartonFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class BudnevEPA : public IntegratedPartonFlux {
  public:
    explicit BudnevEPA(const ParametersList& params) : IntegratedPartonFlux(params) {}

    bool fragmenting() const override final { return false; }
    pdgid_t partonPdgId() const override final { return PDG::photon; }

    static ParametersDescription description() {
      auto desc = IntegratedPartonFlux::description();
      desc.setDescription("Generic Budnev integrated EPA flux");
      desc.add<Limits>("q2range", {0., 1.e4}).setDescription("kinematic range for the parton virtuality, in GeV^2");
      return desc;
    }
  };

  class BudnevEPALepton final : public BudnevEPA {
  public:
    explicit BudnevEPALepton(const ParametersList& params)
        : BudnevEPA(params), ml2_(std::pow(PDG::get().mass(steer<pdgid_t>("pdgId")), 2)) {
      CG_DEBUG("BudnevEPALepton") << "Budnev EPA for photon-from-lepton elastic limit (lepton: "
                                  << PDG::get().name(steer<int>("pdgId")) << ").\n\t "
                                  << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
    }

    static ParametersDescription description() {
      auto desc = BudnevEPA::description();
      desc.setDescription("Lepton integrated EPA (Budnev)");
      desc.add<pdgid_t>("pdgId", PDG::electron).setDescription("lepton PDG id");
      return desc;
    }

    double mass2() const override { return ml2_; }
    double flux(double x) const override {
      Limits q2range;
      if (!computeQ2range(x, q2range))
        return 0.;
      return std::max(0.,
                      prefactor_ * (ml2_ * x * (1. / q2range.max() - 1. / q2range.min()) +
                                    (1. - x + 0.5 * x * x) / x * log(q2range.max() / q2range.min())));
    }

  protected:
    const double ml2_;
  };

  class BudnevEPANucleon : public BudnevEPA {
  public:
    explicit BudnevEPANucleon(const ParametersList& params)
        : BudnevEPA(params), q20_(steer<double>("q20")), mup_(steer<double>("mup")) {}

    static ParametersDescription description() {
      auto desc = BudnevEPA::description();
      desc.add<double>("q20", 0.71);
      desc.add<double>("mup", 2.79).setDescription("proton magnetic moment");
      return desc;
    }

    double flux(double x) const override {
      init();
      Limits q2range;
      if (!computeQ2range(x, q2range))
        return 0.;
      return std::max(0., prefactor_ * (phi(x, q2range.max() / q20_) - phi(x, q2range.min() / q20_)) * (1 - x) / x);
    }

  protected:
    double phi(double x, double qq) const {
      if (x == 1.)
        return 0.;
      const double qq1 = 1 + qq, y = x * x / (1 - x);
      double fac1 = -log(qq1 / qq), fac2 = log((qq1 - b_) / qq1);
      for (size_t i = 1; i <= 3; ++i) {
        const auto fac = std::pow(1. / qq1, i) / i;
        fac1 += fac;
        fac2 += std::pow(b_, i) * fac;
      }
      return (1 + a_ * y) * fac1 + (1 - b_) * y / (4 * qq * qq1 * qq1 * qq1) + c_ * (1 + y / 4) * fac2;
    }

  private:
    void init() const {
      if (init_)
        return;
      const auto mup2 = mup_ * mup_, tau = 4. * mass2() / q20_;
      a_ = 0.25 * (1. + mup2) + tau;
      b_ = 1. - tau;
      c_ = (mup2 - 1.) * std::pow(b_, -4);
      init_ = true;
    }

    const double q20_, mup_;
    mutable bool init_{false};
    mutable double a_{0.}, b_{0.}, c_{0.};
  };

  struct BudnevEPAProton final : public BudnevEPANucleon {
    explicit BudnevEPAProton(const ParametersList& params) : BudnevEPANucleon(params) {
      CG_DEBUG("BudnevEPAProton") << "Budnev EPA for photon-from-proton elastic limit.\n\t"
                                  << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
    }

    static ParametersDescription description() {
      auto desc = BudnevEPANucleon::description();
      desc.setDescription("Proton integrated EPA (Budnev)");
      desc.add<pdgid_t>("pdgId", PDG::proton);
      return desc;
    }
    double mass2() const override { return mp2_; }
  };

  class BudnevEPAHI final : public BudnevEPANucleon {
  public:
    explicit BudnevEPAHI(const ParametersList& params)
        : BudnevEPANucleon(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          mass2_(hi_.mass() * hi_.mass()) {
      CG_DEBUG("BudnevEPAHI") << "Budnev EPA for photon-from-heavy ion elastic limit (HI: " << hi_ << ").\n\t"
                              << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      if (q2_range_.max() < q2max_min_)
        throw CG_FATAL("BudnevEPAHI") << "Invalid maximal Q^2. Its value should not be below " << q2max_min_
                                      << " GeV^2. It is currently " << q2_range_ << ".";
    }

    static ParametersDescription description() {
      auto desc = BudnevEPANucleon::description();
      desc.setDescription("HI integrated EPA (Budnev)");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb()).setDescription("type of heavy ion considered");
      desc.add<Limits>("q2range", {0., 1.e5});
      return desc;
    }

    double mass2() const override { return mass2_; }
    double flux(double x) const override {
      const auto z = (unsigned short)hi_.Z;
      return z * BudnevEPANucleon::flux(x);
    }

  private:
    const HeavyIon hi_;
    const double mass2_, q2max_min_{1.e4};
  };
}  // namespace cepgen
REGISTER_INTEGRATED_PARTON_FLUX("BudnevEPALepton", BudnevEPALepton);
REGISTER_INTEGRATED_PARTON_FLUX("BudnevEPAHI", BudnevEPAHI);
REGISTER_INTEGRATED_PARTON_FLUX("BudnevEPAProton", BudnevEPAProton);
