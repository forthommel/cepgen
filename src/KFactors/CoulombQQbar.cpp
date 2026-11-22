/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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
#include "CepGen/KFactors/KFactor.h"
#include "CepGen/Modules/KFactorFactory.h"
#include "CepGen/Physics/PDG.h"

using namespace std::string_literals;

namespace cepgen::kfactor {
  constexpr double kCF = 4. / 3.;  ///< QCD colour factor of the fundamental representation

  /// Near-threshold (Coulomb/Sommerfeld) QCD correction to the colour-singlet \f$\gamma\gamma\to q\bar q\f$ production
  /// \note At NLO and close to threshold, the leading correction reads
  ///  \f$K=1+\frac{C_F\pi\alpha_S(\mu)}{2\beta}\f$, with \f$\beta=\sqrt{1-4m_q^2/\hat s}\f$ the quark velocity in the pair
  ///  rest frame and \f$C_F=4/3\f$. The \a resum option replaces it with the all-order Sommerfeld factor
  ///  \f$X/(1-e^{-X})\f$, \f$X=C_F\pi\alpha_S/\beta\f$ (whose first-order expansion is the above).
  ///  It is only a threshold approximation, neglecting the (non-singular) constant and high-energy terms of the full
  ///  NLO result: the validity range \a beta is therefore restricted by default. For the complete NLO
  ///  dependence, use the "functional" K-factor with an external reference expression.
  class CoulombQQbar final : public KFactor {
  public:
    explicit CoulombQQbar(const ParametersList& params)
        : KFactor(params), beta_range_(steer<Limits>("beta")), resum_(steer<bool>("resum")) {}

    static ParametersDescription description() {
      auto desc = KFactor::description();
      desc.setDescription("Coulomb threshold correction to gamma gamma -> q qbar");
      desc.add("beta", Limits{0.01, 0.3}).setDescription("validity range on the quark velocity in the pair rest frame");
      desc.add("resum", false).setDescription("resum the Coulomb terms to all orders (Sommerfeld factor)?");
      return desc;
    }

  protected:
    void checkApplicability(const Event& event) const override {
      const auto& central = event(Particle::Role::CentralSystem);
      if (central.size() != 2 || PDG::get().colours(central.at(0).integerPdgId()) <= 1 ||
          central.at(0).integerPdgId() != -central.at(1).integerPdgId())
        throw CG_FATAL("KFactor:coulomb") << "This K-factor requires a central system composed of a quark pair.";
      for (const auto role : {Particle::Role::Parton1, Particle::Role::Parton2})
        if (!event.hasRole(role) || event.oneWithRole(role).pdgId() != PDG::photon)
          throw CG_FATAL("KFactor:coulomb") << "This K-factor requires colour-neutral two-photon initial state.";
    }

    bool inValidityRange(const Event& event) const override { return beta_range_.contains(beta(event)); }

    double compute(const Event& event, double mu) const override {
      const auto b = beta(event);
      const auto x = kCF * M_PI * alphaS(mu) / b;
      return resum_ ? x / (1. - std::exp(-x)) : 1. + 0.5 * x;
    }

  private:
    static double beta(const Event& event) {
      const auto& central = event(Particle::Role::CentralSystem);
      const auto shat = (central.at(0).momentum() + central.at(1).momentum()).mass2();
      return shat > 0. ? std::sqrt(std::max(0., 1. - 4. * central.at(0).momentum().mass2() / shat)) : 0.;
    }

    const Limits beta_range_;
    const bool resum_;
  };
}  // namespace cepgen::kfactor
using cepgen::kfactor::CoulombQQbar;
REGISTER_KFACTOR("coulomb", CoulombQQbar);
