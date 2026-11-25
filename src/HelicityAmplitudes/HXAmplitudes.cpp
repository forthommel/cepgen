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
#include <map>
#include <tuple>

#include "CepGen/Core/Exception.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"

using namespace cepgen;
using namespace cepgen::amplitudes;

/// Model-independent partial-wave helicity amplitudes for \f$\gamma^*\gamma^*\to HX\f$, with \f$H\f$ a spin-0
/// (parity-even) Higgs boson and \f$X\f$ a neutral, C-even particle of integer spin \f$s\f$ and intrinsic parity
/// \f$\eta_X\f$. Only transverse photons are considered.
///
/// Conventions (as in the other models: photon 1 along +z, photon 2 along -z, scattering plane azimuth set to 0):
///  - \f$\Lambda=\lambda_1-\lambda_2\f$ is the spin projection on the beam axis; \f$\lambda_X\f$ the helicity of X,
///    so that the projection on the H direction is \f$-\lambda_X\f$
///  - \f$A_{\lambda_1\lambda_2;\lambda_X}(\theta)=\sum_J(2J+1)a^J_{\lambda_1\lambda_2;\lambda_X}\,d^J_{\Lambda,-\lambda_X}(\theta)\f$
///  - parity: \f$a^J_{-\lambda_1,-\lambda_2;-\lambda_X}=\eta_X(-1)^s a^J_{\lambda_1\lambda_2;\lambda_X}\f$ (J-independent,
///    up to a convention-dependent phase), imposed automatically on the user-provided amplitudes
///  - the amplitudes are returned with a unit coupling: the overall scale is to be applied by the process
class HXAmplitudes final : public HelicityAmplitudes {
public:
  explicit HXAmplitudes(const ParametersList& params)
      : HelicityAmplitudes(params),
        spin_(steer<int>("spin")),
        eta_(steer<int>("parity")),
        jmax_(steer<int>("jmax")),
        landau_yang_(steer<bool>("landauYang")) {
    if (spin_ < 0)
      throw CG_FATAL("HXAmplitudes") << "X must have a non-negative integer spin (single production of a half-integer "
                                     << "spin particle from γγ → H X would violate fermion number conservation).";
    if (eta_ != 1 && eta_ != -1)
      throw CG_FATAL("HXAmplitudes") << "The intrinsic parity of X must be ±1, got " << eta_ << ".";
    if (jmax_ < 0)
      throw CG_FATAL("HXAmplitudes") << "The maximal partial wave must be positive, got " << jmax_ << ".";
    for (const auto& amplitude : steer<std::vector<ParametersList> >("amplitudes")) {
      const auto l1 = amplitude.get<int>("l1"), l2 = amplitude.get<int>("l2");
      if (std::abs(l1) != 1 || std::abs(l2) != 1)
        throw CG_FATAL("HXAmplitudes") << "Only transverse photons are supported (l1, l2 = ±1).";
      const auto lx = amplitude.get<int>("lX");
      if (std::abs(lx) > spin_)
        throw CG_FATAL("HXAmplitudes") << "|λ_X| = " << std::abs(lx) << " exceeds the spin of X (" << spin_ << ").";
      const auto j = amplitude.get<int>("J");
      if (j < 0 || j > jmax_)
        throw CG_FATAL("HXAmplitudes") << "Partial wave J=" << j << " out of range [0," << jmax_ << "].";
      if (j < std::abs(l1 - l2) || j < std::abs(lx))
        throw CG_FATAL("HXAmplitudes") << "J=" << j << " < max(|Λ|,|λ_X|) is not allowed.";
      if (landau_yang_ && l1 == l2 && j % 2 == 1)
        throw CG_FATAL("HXAmplitudes") << "Landau-Yang: odd J is forbidden for λ₁=λ₂.";
      amplitudes_[{j, +l1, +l2, +lx}] = amplitude.get<std::complex<double> >("amplitude", 0.);
      amplitudes_[{j, -l1, -l2, -lx}] =
          amplitudes_.at({j, +l1, +l2, +lx}) * (eta_ * (spin_ % 2 ? -1. : 1.));  // parity partner
    }
  }

  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Model-independent γγ → H X partial-wave helicity amplitudes");
    desc.add("spin", 0).setDescription("integer spin of X");
    desc.add("parity", +1).setDescription("intrinsic parity of X (±1); natural parity if equal to (-1)^spin");
    desc.add("jmax", 0).setDescription("highest partial wave included");
    desc.add("landauYang", false).setDescription("forbid odd J for λ₁=λ₂ (strictly valid for on-shell photons only)");
    desc.add("amplitudes", std::vector<ParametersList>{})
        .setDescription("partial-wave amplitudes {J, l1, l2, lX, re, im}; parity partners are implicit");
    return desc;
  }

  bool hasAAHX() const override { return true; }

  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const int l1 = hel[0], l2 = hel[1], lh = hel[2], lx = hel[3];
    if (std::abs(l1) != 1 || std::abs(l2) != 1 || lh != 0 || std::abs(lx) > spin_)  // transverse photons, spin-0 H
      return 0.;
    const auto theta = std::acos(kin.cos_theta);
    const int lambda = l1 - l2;
    std::complex<double> out = 0.;
    for (int j = 0; j <= jmax_; ++j)
      if (const auto it = amplitudes_.find({j, l1, l2, lx}); it != amplitudes_.end())
        out += (2. * j + 1.) * it->second * wignerD(j, lambda, -lx, theta);
    return out;
  }

private:
  /// Wigner small-d function \f$d^j_{m'm}(\beta)\f$ for integer arguments
  static double wignerD(int j, int mp, int m, double beta) {
    if (std::abs(mp) > j || std::abs(m) > j)
      return 0.;
    const auto fact = [](int n) { return std::tgamma(n + 1.); };
    const auto c = std::cos(0.5 * beta), s = std::sin(0.5 * beta);
    const auto pref = std::sqrt(fact(j + mp) * fact(j - mp) * fact(j + m) * fact(j - m));
    double sum = 0.;
    for (int k = 0; k <= 2 * j + 1; ++k) {
      const int a = j + m - k, b = mp - m + k, d = j - mp - k;
      if (a < 0 || b < 0 || d < 0)
        continue;
      sum += (((mp - m + k) % 2 == 0) ? 1. : -1.) * std::pow(c, 2 * j + m - mp - 2 * k) * std::pow(s, mp - m + 2 * k) /
             (fact(a) * fact(k) * fact(b) * fact(d));
    }
    return pref * sum;
  }

  const int spin_, eta_, jmax_;
  const bool landau_yang_;
  std::map<std::tuple<int, int, int, int>, std::complex<double> > amplitudes_;
};
REGISTER_HELICITY_AMPLITUDES("hxModel", HXAmplitudes);
