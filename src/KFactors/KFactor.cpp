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
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

using namespace cepgen;
using namespace std::string_literals;

namespace {
  auto parseScale(const std::string& name) {
    using Scale = KFactor::Scale;
    if (name == "fixed")
      return Scale::fixed;
    if (name == "mass")
      return Scale::mass;
    if (name == "pt")
      return Scale::pt;
    throw CG_FATAL("KFactor") << "Invalid scale choice: '" << name << "'. Allowed: fixed, mass, pt.";
  }
}  // namespace

KFactor::KFactor(const ParametersList& params)
    : NamedModule(params),
      scale_(parseScale(steer<std::string>("scale"))),
      mu_(steer<double>("mu")),
      mu_factor_(steer<double>("muFactor")),
      mu_min_(steer<double>("muMin")),
      mass_range_(steer<Limits>("mass")),
      reject_(steer<std::string>("outOfRange") == "reject"),
      alpha_s_(AlphaSFactory::get().build(steer<ParametersList>("alphaS"))) {
  if (const auto policy = steer<std::string>("outOfRange"); policy != "unity" && policy != "reject")
    throw CG_FATAL("KFactor") << "Invalid out-of-range policy: '" << policy << "'. Allowed: unity, reject.";
  if (mu_factor_ <= 0. || mu_min_ <= 0.)
    throw CG_FATAL("KFactor") << "muFactor and muMin must be strictly positive.";
  if (scale_ == Scale::fixed && mu_ <= 0.)
    throw CG_FATAL("KFactor") << "A strictly positive fixed scale 'mu' is required with scale='fixed'.";
}

KFactor::~KFactor() = default;

ParametersDescription KFactor::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed higher-order correction (K-factor) modelling");
  desc.add("scale", "mass"s).setDescription("renormalisation scale base choice (fixed, mass, pt)");
  desc.add("mu", 0.).setDescription("fixed renormalisation scale (GeV), if scale=fixed");
  desc.add("muFactor", 1.).setDescription("multiplicative factor on the scale (e.g. 0.5, 2 for scale variations)");
  desc.add("muMin", 1.).setDescription("lowest renormalisation scale allowed (GeV)");
  desc.add("mass", Limits{}).setDescription("validity range on the central system invariant mass (GeV)");
  desc.add("outOfRange", "unity"s).setDescription("action outside validity range (unity: K=1, reject: K=0)");
  desc.add("alphaS", AlphaSFactory::get().describeParameters("pegasus"))
      .setDescription("strong coupling evolution algorithm");
  return desc;
}

bool KFactor::singleThreaded() const { return alpha_s_ && alpha_s_->singleThreaded(); }

double KFactor::alphaS(double mu) const { return (*alpha_s_)(mu); }

double KFactor::scale(const Event& event) const {
  double mu0 = mu_;
  if (scale_ != Scale::fixed) {
    const auto& central = event(Particle::Role::CentralSystem);
    if (scale_ == Scale::mass) {
      Momentum sum;
      for (const auto& part : central)
        sum += part.momentum();
      mu0 = sum.mass();
    } else {
      mu0 = 0.;
      for (const auto& part : central)
        mu0 += part.momentum().pt();
      mu0 /= std::max<size_t>(1, central.size());
    }
  }
  return std::max(mu_factor_ * mu0, mu_min_);
}

double KFactor::operator()(const Event& event) const {
  if (!checked_.exchange(true))
    checkApplicability(event);
  Momentum sum;
  for (const auto& part : event(Particle::Role::CentralSystem))
    sum += part.momentum();
  if (!mass_range_.contains(sum.mass()) || !inValidityRange(event)) {
    if (!warned_.exchange(true))
      CG_WARNING("KFactor") << "K-factor '" << name() << "' evaluated outside of its validity range for the first "
                            << "time (central mass: " << sum.mass() << " GeV). Applying the '"
                            << (reject_ ? "reject" : "unity") << "' policy for such points.";
    return reject_ ? 0. : 1.;
  }
  return compute(event, scale(event));
}
