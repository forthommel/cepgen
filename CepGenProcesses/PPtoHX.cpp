/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

/// Define X in the PDG library (before the base class is built) and return the central system content
static spdgids_t centralIds(const ParametersList& params) {
  const auto id = params.get<int>("pdgIdX", 9900025);
  PDG::get().define(
      ParticleProperties(static_cast<pdgid_t>(id), "X", "X", 0, params.get<double>("massX", 100.), 0., {0}, false));
  return spdgids_t{static_cast<spdgid_t>(PDG::higgs), static_cast<spdgid_t>(id)};
}

/// Compute the matrix element for a CE \f$\gamma^*\gamma^*\rightarrow H X\f$ process, with \f$X\f$ a neutral, stable particle
/// \note The spin, parity and \f$\gamma\gamma HX\f$ couplings of X are described by a helicity amplitudes model
///  (default: "hxModel", model-independent partial waves)
class PPtoHX final : public proc::FactorisedProcess {
public:
  explicit PPtoHX(const ParametersList& params)
      : proc::FactorisedProcess(params, centralIds(params)), mx_(steer<double>("massX")) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoHX>(parameters()); }

  static ParametersDescription description() {
    auto desc = proc::FactorisedProcess::description();
    desc.setDescription("γ*γ* → H X (X neutral and stable, helicity amplitudes model)");
    desc.add("ktFactorised", true);
    desc.add("randomiseCharge", false);  // H and X are neutral: no charge conjugation to randomise
    desc.add("pdgIdX", 9900025).setDescription("PDG id assigned to X (defined on the fly in the PDG library)");
    desc.add("massX", 100.).setDescription("mass of X (GeV/c²)");
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("hxModel"))
        .setDescription("Helicity amplitudes model");
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    if (!partonsKtFactorised())
      throw CG_FATAL("PPtoHX:prepare") << "This process is only defined for kt-factorised (off-shell) photons.";
    if (mh2_ = std::pow(PDG::get()(PDG::higgs).mass, 2); mh2_ <= 0.)
      throw CG_FATAL("PPtoHX:prepare") << "The Higgs mass is not defined in the PDG library.";
    mx2_ = mx_ * mx_;
    if (amplitudes_ = amplitudes::HelicityAmplitudesFactory::get().build(steer<ParametersList>("model"));
        !amplitudes_->hasAAHX())
      throw CG_FATAL("PPtoHX:prepare") << "Model '" << amplitudes_->name()
                                       << "' does not have a registered γγ → H X vertex in its specifications.";
    spin_ = amplitudes_->parameters().get<int>("spin", 0);
    CG_INFO("PPtoHX:prepare") << "γγ → H X: mH = " << std::sqrt(mh2_) << " GeV, mX = " << mx_ << " GeV. "
                              << "Helicity amplitudes (model: " << amplitudes_->name()
                              << ") initialised. Parameters: " << amplitudes_->parameters() << ".";
  }

  double computeFactorisedMatrixElement() override {
    if (q1().pt2() == 0. || q2().pt2() == 0.)
      return 0.;
    const auto s_hat = shat();
    if (s_hat <= std::pow(std::sqrt(mh2_) + mx_, 2))
      return 0.;
    const amplitudes::Kinematics kin(mh2_, mx2_, s_hat, that(), uhat());
    // photon polarisations are taken relative to the central pair azimuth (amplitudes defined in the scattering plane)
    const auto phi_pair = std::atan2(pc(0).py() - pc(1).py(), pc(0).px() - pc(1).px());
    const auto e1 = std::polar(1., std::atan2(q1().py(), q1().px()) - phi_pair),
               e2 = std::polar(1., std::atan2(q2().py(), q2().px()) - phi_pair);
    double amat2{0.};
    for (short lam_x = -spin_; lam_x <= spin_; ++lam_x)  // H has helicity 0
      amat2 += std::norm(amplitudes_->amplitudeLinearPolarisations(kin, 0, lam_x, e1, e2));
    if (!utils::positive(amat2))
      return 0.;
    // the model returns the full dimensionless amplitude M (contact terms and resonances): |M|²/ŝ² matches the
    // 2→2 phase space convention
    return amat2 / (s_hat * s_hat);
  }

  const double mx_;

  double mh2_{0.};
  double mx2_{0.};
  short spin_{0};
  std::unique_ptr<const amplitudes::HelicityAmplitudes> amplitudes_;
};
REGISTER_PROCESS("pptohx", PPtoHX);
